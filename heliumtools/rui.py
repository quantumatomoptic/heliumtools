from heliumtools.correlations import Correlation
from heliumtools.tools import apply_ROI
from heliumtools.misc.some_plots_volume1 import heatmap_with_boxes

from scipy.optimize import curve_fit
from scipy.interpolate import RegularGridInterpolator
from scipy.special import erf
from scipy.integrate import quad

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import copy

class Correlation1D(Correlation):

    def __init__(self, atoms, **kwargs):
        super().__init__(atoms, **kwargs)
        self.var1 = None
        self.var2 = None
        self.round_decimal = 7
        self.boxes = {
            "1": {
                "Vx": {"size": 10, "position": 0},
                "Vy": {"size": 3, "position": 0},
                "Vz": {"size": 0.9, "position": 56},
            },
            "2": {
                "Vx": {"size": 10, "position": 0},
                "Vy": {"size": 3, "position": 0},
                "Vz": {"size": 0.9, "position": 130},
            },
        }
        self.is_there_a_copy_of_total = False
        self.compute_errors = False
        self.remove_shot_noise = True
        self.__dict__.update(kwargs)
        self.boxes = copy.deepcopy(self.boxes)
        self.width = [11, 14]
        self.length = 3.0
        self.sliding = 3
        self.step_size = 0.5
        self.correlation_order_max = 4
        self.regions = ["cross","loc1","loc2"]

    def get_atoms_in_boxes_one_variable(self, df_atoms, var, box):
        """ returns a dataframe containing the label of the atoms in the box. 

        Parameters
        ----------
        df_atoms : pd.DataFrame
            atoms dataframe for which the tranvserse selection has already been done.
        var : Variable
            the variable of the box that is scanned, can be box 1 or box 2.
        box : dict
            the box with the axis over which we will iterate.
        
        return
        --------
        atoms_along_ax : pd.DataFrame
            pd.Dataframe with collumns {"Cycle" , "label" , "axis of scan"}. "label" is the label of the atom in a box along the "axis of scan".
            "axis of scan" are the boxes along the scanned axis.
        """

        # 
        atoms_along_ax = []

        # iterate over the box for different values
        for i in range(var.n_step):
            # create box
            box[var.axe][var.type] = var.values[i]
            # get atoms in box
            df = apply_ROI(df_atoms, box)[["label", "Cycle"]]
            # add box to dataframe
            df[var.name] = var.values[i]
            # add dataframe
            atoms_along_ax.append(df)
    
        return pd.concat(atoms_along_ax)
           
    def switch_box_num(self, name):
        """swith the number of box 1 and box 2 in a string. 
        Exemple: switch_box_num('Vz1_Vz2') returns 'Vz2_Vz1'

        Parameters
        ----------
        name : string
            the name for which we wich to change the box number

        Returns
        -------
        string
            the name with number switched
        """
        return name.replace(self.var1.box, "%ù€").replace(self.var2.box,self.var1.box).replace("%ù€", self.var2.box)
    
    def superimpose(self):
        """ This function labels each atom and atributes it to a given box along the scanned axis. It does this in three different regions: local region in the positive momentum 
        region (k,k), local region in the negative momentum region (-k,-k) and finally in the cross region (k,-k). 
        Here k in the momentum that we scan. """

        # middle of the region of integration
        middle = np.mean(self.width)
        # length of integration
        dW =  (np.max(self.width) - np.min(self.width))
        maxi = self.length + dW
        
        # Define boxes to compute 2D map
        self.define_variable1(box="1",axe="Vz",type="position",name="Vz1",min = -middle - maxi,max = -middle + maxi + self.step_size,step = self.step_size)
        self.define_variable2(box="2",axe="Vz",type="position",name="Vz2",min = middle - maxi,max = middle + maxi + self.step_size,step = self.step_size)

        # first we label all atoms
        self.atoms.reset_index(drop=True)
        self.atoms["label"] = np.arange(0, len(self.atoms))

        """  Get atoms in the first box var1 """
        # create copy of box defined in var1
        box = self.boxes[self.var1.box].copy() # transverse selection
        # remove the scanned axis of the box, e.g Vz
        posi_and_size = box.pop(self.var1.axe)
        # create scanned box
        scanned_box = {self.var1.axe: posi_and_size}
        # get all atoms that are inside the transverse box
        df_atoms = apply_ROI(self.atoms, box)
        # seperate the atoms into boxes along the axis of scan
        self.df_atoms_var1 = self.get_atoms_in_boxes_one_variable(df_atoms, self.var1, scanned_box)
        # round value of box
        self.df_atoms_var1[self.var1.name] = np.round(self.df_atoms_var1[self.var1.name],4)
        
        """  Get atoms in the second box var2. Same as for var1 """
        box = self.boxes[self.var2.box].copy() # transverse selection
        posi_and_size = box.pop(self.var1.axe)
        scanned_box = {self.var2.axe: posi_and_size}
        df_atoms = apply_ROI(self.atoms, box)
        self.df_atoms_var2 = self.get_atoms_in_boxes_one_variable(df_atoms, self.var2, scanned_box,)
        self.df_atoms_var2[self.var2.name] = np.round(self.df_atoms_var2[self.var2.name],4)

        """ reset index """
        self.df_atoms_var1.reset_index(drop=True,inplace=True)
        self.df_atoms_var2.reset_index(drop=True,inplace=True)

    def compute_correlations_superimpose(self):
        """ From the dataframes built in self.compute_superimpose(), this function calculates correlations in three different regions: local correlation
        in the positive momentum region (k,k), local correlation in the negative momentum region (-k,-k) and finally cross corrlation in the cross region (k,-k). 
        Here k is the momentum that we scan. In the end for each region we have two kinds of dataset: one where we have the atoms in each box in each cycle
        and another where we have number of atoms per box """

        """ compute atoms correlation dataframe where each atom is assigned to a box """
        
        # merge datasets to compute correlation in cross region (-k,k)
        self.cross = self.df_atoms_var2.merge(self.df_atoms_var1, how="outer", on="Cycle").dropna()
        # merge datasets to compute correlation in cross region (-k,-k)
        self.loc1 = self.df_atoms_var1.merge(self.df_atoms_var1.rename(columns = {self.var1.name:self.var2.name}), how="outer", on="Cycle").dropna()
        # merge datasets to compute correlation in cross region (k,k)
        self.loc2 = self.df_atoms_var2.merge(self.df_atoms_var2.rename(columns = {self.var2.name:self.var1.name}), how="outer", on="Cycle").dropna()

        """ compute atom number in each box (summed over cycles) """

        # count atom number in each box along the scanned axis
        self.density_var1 = self.df_atoms_var1.groupby(self.var1.name).count().reset_index()[[self.var1.name, "label"]].rename(columns={'label': "N_" + self.var1.box})
        self.density_var2 = self.df_atoms_var2.groupby(self.var2.name).count().reset_index()[[self.var2.name, "label"]].rename(columns={'label': "N_" + self.var2.box})

        """ Merge dataframes self.density_var1 and self.density_var2 in the three different regions: (-k,k), (-k,-k) and (k,k) """

        # merge dataframes to compute density in cross region (-k,k)
        self.cross_result = pd.merge(self.density_var2, self.density_var1, how = "cross")
        # merge dataframes to compute density in local region (-k,-k)
        self.loc1_result = pd.merge(self.density_var1,self.density_var1.rename(columns={"N_" + self.var1.box:"N_" + self.var2.box,self.var1.name:self.var2.name}),how = "cross")
        # merge dataframes to compute density in local region (k,k)
        self.loc2_result = pd.merge(self.density_var2,self.density_var2.rename(columns={"N_" + self.var2.box:"N_" + self.var1.box,self.var2.name:self.var1.name}),how = "cross")

        """ 2D correlation map """

        # we start by counting the atoms in the 2D boxes [var1,var2] for each region
        for attr_result, attr_atoms in zip(["cross_result", "loc1_result", "loc2_result"],["cross", "loc1", "loc2"]):
            # get density dataframe ["cross_result", "loc1_result", "loc2_result"] from the object.
            df_result = getattr(self, attr_result)
            # get atom dataframe ["cross", "loc1", "loc2"] from the object.
            df_atoms = getattr(self, attr_atoms)

            # we are carefull to not count the same atom twice, by only selecting rows with atoms different labels corresponding to different axis var1 and Var2.
            # we count the atom number in each 2D box with sides [var1,var2]
            df_at = (
                df_atoms[df_atoms['label_x'] != df_atoms['label_y']]
                .groupby([self.var2.name, self.var1.name])
                .count()
                .reset_index()
                .rename(columns={"Cycle": "G2"})
                .drop(columns=["label_x", "label_y"])
            )

            # merge dataframes. We now have in the same dataframe, the density and correlations
            df_result =  df_result.merge(df_at, on=[self.var2.name, self.var1.name], how="outer")

            # compute mean atom number in boxes
            df_result["N_" + self.var2.box] = df_result["N_" + self.var2.box]/self.n_cycles
            df_result["N_" + self.var1.box] = df_result["N_" + self.var1.box]/self.n_cycles
            
            # calculate G2 and g2
            df_result["G2"] = df_result["G2"]/self.n_cycles
            df_result["g^2"] = df_result["G2"]/ (df_result["N_" + self.var2.box]*df_result["N_" + self.var1.box])

            # replace NaN by zeros
            df_result = df_result.replace(np.nan, 0.0)

            df_result["<N1><N2>"] = df_result["N_1"]*df_result["N_2"]

            # Assign the merged DataFrame back to the attribute
            setattr(self, attr_result,df_result)
        
    def get_error(self,df, columns = ["N_1","N_2"], stat = "thermal", max_val = 100):
        """return the relative error assuming thermal or poissonian statistics depending on the number of variables

        Parameters
        ----------
        df : pandas.DataFrame
            the datarame that contains the data
        columns : list, optional
            the keys on which the error mus be computed, by default ["N_1","N_2"]
        stat : str, optional
            the statistics of the error, thermal or poissonian, by default "thermal"
        max_val : int, optional
            the value to put if error is Nan, by default 100

        Returns
        -------
        pandas.Series
            the relative error
        """
        err = np.zeros(len(df))
        for col in columns:
            if stat in "thermal":
                err += (df[col]**2+df[col])/df[col]**2
            else:## poisson
                err += 1/df[col]
        err = np.nan_to_num(np.sqrt(err)/np.sqrt(self.n_cycles),nan=max_val, posinf=max_val, neginf=max_val)
        return err
    
    def compute_1D_correlation(self):
        """ From the dataframe with each atom in a box, we compute the 1D correlation by integrating the 2D map in a region along the diagonal.
        Parameters
        --------------
        self.width : array of size 2
            lower and upper limit of integration
        self.length: float
            half the size to compute 1D correlation along diagonal: 1D correlation from -corr.length to corr.length
        """
        
        # middle of the region of integration
        middle = np.mean(self.width)
        # length of integration
        dW =  (np.max(self.width) - np.min(self.width))
        maxi = self.length + dW

        # define diagonal and anti-diagonals
        roi_list = [{"Vz1-Vz2": {"position": -2 * middle, "size": 2*dW}},
                    {"Vz1+Vz2": {"position": -2 * middle, "size": 2*dW}},
                    {"Vz1+Vz2": {"position": 2 * middle, "size": 2*dW}}]
           
        # compute 1D correlation along diagonal and anti diagonals
        for j, (attr_result, df_atoms, df_initial_res,axis_key) in enumerate(zip(
            ["cross_1Dresult", "loc1_1Dresult", "loc2_1Dresult"],
            [self.cross, self.loc1,self.loc2], 
            [self.cross_result, self.loc1_result, self.loc2_result],
            ["Vz1+Vz2","Vz1-Vz2","Vz1-Vz2"],
            )):

            """ We compute diagonals and round the values to avoid numerical problems """
            # on dataframe with atom number per box
            df_initial_res["Vz1+Vz2"]=round(df_initial_res["Vz1"]+df_initial_res["Vz2"],3)
            df_initial_res["Vz1-Vz2"]=round(df_initial_res["Vz1"]-df_initial_res["Vz2"],3)
            # on dataframe with atoms
            df_atoms["Vz1+Vz2"] = round(df_atoms["Vz1"]+df_atoms["Vz2"],3)
            df_atoms["Vz1-Vz2"] = round(df_atoms["Vz1"]-df_atoms["Vz2"],3)

            # round values to avoid numerical problems
            df_initial_res["Vz1"] = np.round(df_initial_res["Vz1"],3)
            df_initial_res["Vz2"] = np.round(df_initial_res["Vz2"],3)
            df_initial_res["<N_1><N_2>"] = df_initial_res["N_1"]*df_initial_res["N_2"]
                
            """ Apply the selection along the integration axis """
            # on dataframe with atom number per box
            df_initial_res_roi = apply_ROI(df_initial_res, roi_list[j])
            # on dataframe with atoms
            df_atoms_roi = apply_ROI(df_atoms, roi_list[j])

            """ Integration over the integrated axis """
            
            # on dataframe with atom number per box
            df_1Dresult = df_initial_res_roi.groupby(axis_key).mean().reset_index() # this df contains g^2
            df_1Dsum = df_initial_res_roi.groupby(axis_key).sum().reset_index()
            df_1Dresult["g^(2) integrated"] = df_1Dresult["G2"] / df_1Dresult["<N_1><N_2>"]

            # on dataframe with atom labels
            df_at = (
                df_atoms_roi[df_atoms_roi['label_x'] != df_atoms_roi['label_y']]
                .groupby(axis_key)
                .count()
                .reset_index()
                .rename(columns={"Cycle": "Int(G2)"})[[axis_key,"Int(G2)" ]]
            )
            df_at["Int(G2)"] = df_at["Int(G2)"]/self.n_cycles

            """ merge dataframes """
            df_1Dresult = df_1Dresult.merge(df_at, on = axis_key, how = "outer")
            df_1Dresult = df_1Dresult.merge(df_1Dsum.rename(columns={"N_1": "Int(N_1)",
                                                                     "<N_1><N_2>": "Int(<N_1><N_2>)",
                                                                     "N_2": "Int(N_2)"},
                                                                     )[["Int(N_1)", "Int(N_2)", axis_key,"Int(<N_1><N_2>)"]],
                                                                         on = axis_key, how="outer")
            df_1Dresult["g^(2) int2"] = df_1Dresult["Int(G2)"]/df_1Dresult["Int(<N_1><N_2>)"]

            ## Assign the result to self (the class)
            ## so that the name of the final dataframe will be 
            # self.cross_1Dresult, self.loc1_1Dresult and  self.loc2_1Dresult
            setattr(self, attr_result,df_1Dresult)
            
    def fit1DCorrelations(self):
        """ Function that fits local and cross normalized correlations 
        Return
        --------------
        fitRes : dictionary of dictionaries
            dataset with fit results
        """

        # create dataset to hold fit results
        fitRes = dict()
        keys = ["cross","loc1","loc2"]
        for key in keys:
            fitRes[key] = dict()

        # for each region fit correlation
        for j, (attr_result,axis_key,key) in enumerate(zip([self.cross_1Dresult, self.loc1_1Dresult, self.loc2_1Dresult],["Vz1+Vz2","Vz1-Vz2","Vz1-Vz2"],keys)):
            # get data in region of interest
            my_df = apply_ROI(attr_result, {axis_key: [-self.length, self.length]}).dropna()
            # for each way of calculating g2
            for jjj, Ystr in enumerate(["g^2", "g^(2) integrated","g^(2) int2"]):
                # calculate error of measurement
                x = my_df[axis_key]
                y = my_df[Ystr]
                err = np.nan_to_num(np.sqrt(my_df["G2"]/my_df["G2"]**2+my_df["N_1"]/my_df["N_1"]**2+my_df["N_2"]/my_df["N_2"]**2)/np.sqrt(self.n_cycles), nan=100, posinf=100, neginf=100)
                # try to fit gaussian
                try: 
                    fit_func = gaussian
                    popt, pcov = curve_fit(fit_func, x, y, p0=[1, 1, 0], sigma=err,absolute_sigma=True, bounds = ([0,0,-0.2],[5,3,0.2]))
                    maxi  = popt[0]+1
                    perr = np.sqrt(np.diag(pcov))[0]

                    # add fit resuts to dataset
                    fitRes[key][Ystr] = dict()
                    fitRes[key][Ystr]["value"] = maxi
                    fitRes[key][Ystr]["error"] = perr
                except Exception as e:
                    print(f"Failed to fit {e}")
                    fitRes[key][Ystr] = dict()
                    fitRes[key][Ystr]["value"] = 0.0
                    fitRes[key][Ystr]["error"] = 0.0
        
        return fitRes
    
    def showCorrelations(self):
        """ Plots the 2D and 1D integrated cross and local correlations in the 3 regions of interest """
        
        # middle of the region of integration
        middle = np.mean(self.width)
        # length of integration
        dW =  (np.max(self.width) - np.min(self.width))
        maxi = self.length + dW
        
        # Define boxes of 2D map
        self.define_variable1(box="1",axe="Vz",type="position",name="Vz1",min = -middle - maxi,max = -middle + maxi + self.step_size,step = self.step_size)
        self.define_variable2(box="2",axe="Vz",type="position",name="Vz2",min = middle - maxi,max = middle + maxi + self.step_size,step = self.step_size)

        # define diagonal and anti-diagonals
        roi_list = [{"Vz1-Vz2": {"position": -2 * middle, "size": 2*dW}},
                    {"Vz1+Vz2": {"position": -2 * middle, "size": 2*dW}},
                    {"Vz1+Vz2": {"position": 2 * middle, "size": 2*dW}}]
        
        # colors and markers for plot
        cmaps = ["Greens", "Blues", "Reds"]
        markers = [ 'o',"v", "s", "d", "p", "*", "H", "P", "<", "+" ]*3

        # initialize plot
        fig, axes = plt.subplots(figsize = (15,15),ncols = 3, nrows=4)
           
        # for each region
        for j, (attr_result, df_atoms, df_initial_res,axis_key, axis_integrated) in enumerate(zip(
            [self.cross_1Dresult, self.loc1_1Dresult, self.loc2_1Dresult],
            [self.cross, self.loc1,self.loc2], 
            [self.cross_result, self.loc1_result, self.loc2_result],
            ["Vz1+Vz2","Vz1-Vz2","Vz1-Vz2"],
            ["Vz1-Vz2","Vz1+Vz2","Vz1+Vz2"],
            )):
            
            # check up what you do, by plotting <N1><N2>, G² and g²
            heatmap_with_boxes(ax=axes[0,j],df=df_initial_res,columns=self.var1.name,index=self.var2.name,boxes={},values="<N1><N2>",cmap="Greys")
            heatmap_with_boxes(ax=axes[1,j],df=df_initial_res,columns=self.var1.name,index=self.var2.name,boxes={},values="G2",cmap="Greys")
            heatmap_with_boxes(ax=axes[2,j],df=df_initial_res,columns=self.var1.name,index=self.var2.name,boxes={},values="g^2",cmap="Greys",vmax=2.2,vmin=0.9)
                
            """ Apply the selection along the integration axis """
            # on dataframe with atom number per box
            df_initial_res_roi = apply_ROI(df_initial_res, roi_list[j])

            # check up what you do, by plotting regions of integration of the 2D plots
            heatmap_with_boxes(ax=axes[0,j],df=df_initial_res_roi,columns=self.var1.name,index=self.var2.name,boxes={},values="<N1><N2>",cmap=cmaps[j],cbar_bool=False,)
            heatmap_with_boxes(ax=axes[1,j],df=df_initial_res_roi,columns=self.var1.name,index=self.var2.name,boxes={},values="G2",cmap=cmaps[j],cbar_bool=False,)
            heatmap_with_boxes(ax=axes[2,j],df=df_initial_res_roi,columns=self.var1.name,index=self.var2.name,boxes={},values="g^2",cmap = cmaps[j],vmax=2.2,vmin=0.9,cbar_bool=False,)
            a = df_initial_res_roi[df_initial_res_roi[axis_integrated] == df_initial_res_roi[axis_integrated].min()]
            b = df_initial_res_roi[df_initial_res_roi[axis_integrated] == df_initial_res_roi[axis_integrated].max()]
            for k in range(3):
                axes[k,j].plot(a["Vz1"],a["Vz2"],color = plt.get_cmap(cmaps[j])(1.0))
                axes[k,j].plot(b["Vz1"],b["Vz2"],color = plt.get_cmap(cmaps[j])(1.0))

           
            # check up what you do, by plotting 1D correlation
            # get rows inside the region of interest
            my_df = apply_ROI(attr_result, {axis_key: [-self.length, self.length]}).dropna()

            ax = axes[3,j]
            for jjj, Ystr in enumerate(["g^2"]):
                color = plt.get_cmap(cmaps[j])(1.4)
                x, y = my_df[axis_key], my_df[Ystr]
                err = np.nan_to_num(np.sqrt(my_df["G2"]/my_df["G2"]**2+my_df["N_1"]/my_df["N_1"]**2+my_df["N_2"]/my_df["N_2"]**2)/np.sqrt(self.n_cycles), nan=100, posinf=100, neginf=100)
                ax.errorbar(x, y,  yerr=err,fmt= markers[jjj], markerfacecolor="none",
                            markeredgecolor =color, ecolor = color)
                xth = np.linspace(np.min(x), np.max(x), 100)
                # try to fit gaussian
                try: 
                    fit_func = gaussian
                    popt, pcov = curve_fit(fit_func, x, y, p0=[1, 1, 0], sigma=err,absolute_sigma=True, bounds = ([0,0,-0.2],[5,3,0.2]))
                    maxi  = popt[0]+1
                    perr = np.sqrt(np.diag(pcov))[0]
                    legenda = "A0 = "+"{:.2f}({:.0f})".format(maxi,100*perr)+"\n max = "+str(np.amax(y))
                    ax.plot(xth, fit_func(xth, *popt), color=color, label = legenda)
                except Exception as e:
                    print(f"Failed to fit {e}")
                ax.grid(True, alpha = 0.5)
                ax.legend()
                ax.set_ylim([.8,2.8])
                ax.set_xlabel(axis_key)

        plt.show()
                
    def fitProdAndG2(self,funcG2,funcProd,guessG2,guessProd,limits,show = False):
        """ Fits a function funcG2 to function G2 and funcProd and <N1><N2> at the same time 
        to ensure that they have the same offset,  in each region.
        Parameters
        --------------
        funcG2 : dictionary of functions used to fit G2
            each function in funcG2 should be a function of two variables, x and y, that should be passed as an array [x,y]. 
            Rest of the arguments for funcG2 are considered to be fit parameters. E.g: funcG2 can be a Skew2DGaussian
        funcProd : disctionary of functions used to fit <N1><N2>
            Same type of as funcG2. E.g: func can be a ProdGaussian
        guessG2 : dictionary of 1D numpy array
            each array is the initial guess for fit parameters of G2, offset included
        guessProd : dictionary of 1D numpy array
            each array is the initial guess for fit parameters of <N1><N2>, offset included
        limits : dictionary of 2D numpy array
            each element of limits is of the form [limit_x,limit_y]. limit_x are the limits of the region of interest of the fit on the x-axis. 
            Equivalent for limit_y on the y-axis.
        show : bool
            True if we want to plot fit result. False otherwise
        """

        # create datastructure to store data
        self.poptG2 = dict()
        self.poptProd = dict()
        for region in self.regions:
            self.poptG2[region] = []
            self.poptProd[region] = []

        # for each region, fit G2 and <N1><N2>
        for j, (region,df) in enumerate(zip(self.regions,
                                            [self.cross_result, self.loc1_result,self.loc2_result])):

            # combined fit of G2 and <N1><N2> so that they have equal offset
            poptG2 , poptProd , pcov = FitG2andProd(df,self.var1.name,self.var2.name,
                                                    funcG2[region],funcProd[region],
                                                    guessG2[region],guessProd[region],limits[region],region,show)
            # store values of fit parameters
            self.poptG2[region] = poptG2
            self.poptProd[region] = poptProd

        # add fit functions to class
        self.funcG2 = funcG2
        self.funcProd = funcProd

    def fitDCEpeaksAndG2(self,funcG2,funcN1,funcN2,guessG2,guessN1,guessN2,limits,show = False):
        """ Fits a function funcG2 to function G2, funcN1 to <N1> and funcN2 to <N2> at the same time 
        to ensure that G2 and <N1><N2> have the same offset,  in each region.
        Parameters
        --------------
        funcG2 : dictionary of functions used to fit G2
            each function in funcG2 should be a function of two variables, x and y, that should be passed as an array [x,y]. 
            Rest of the arguments for funcG2 are considered to be fit parameters. E.g: funcG2 can be a Skew2DGaussian
        funcN1 : dictionary of functions used to fit <N1>
            Same type of as funcG2. E.g: func can be a DCE_N1
        funcN2 : dictionary of functions used to fit <N2>
            Same type of as funcG2. E.g: func can be a DCE_N2
        guessG2 : dictionary of 1D numpy array
            each array is the initial guess for fit parameters of G2, offset included
        guessN1 : dictionary of 1D numpy array
            each array is the initial guess for fit parameters of <N1>, offset included
        guessN2 : dictionary of 1D numpy array
            eacDCEh array is the initial guess for fit parameters of <N2>, offset included
        limits : dictionary of 2D numpy array
            each element of limits is of the form [limit_x,limit_y]. limit_x are the limits of the region of interest of the fit on the x-axis. 
            Equivalent for limit_y on the y-axis.
        show : bool
            True if we want to plot fit result. False otherwise
        """

        # create datastructure to store data
        self.poptG2 = dict()
        self.poptN1 = dict()
        self.poptN2 = dict()
        for region in self.regions:
            self.poptG2[region] = []
            self.poptN1[region] = []
            self.poptN2[region] = []

        # for each region, fit G2 and <N1><N2>
        for j, (region,df) in enumerate(zip(self.regions,
                                            [self.cross_result, self.loc1_result,self.loc2_result])):

            # combined fit of G2, <N1> and <N2> so that G2 and <N1><N2> have equal offset            
            poptG2 , poptN1 , poptN2 , pcov = FitG2andDCEPeaks(df,self.var1.name,self.var2.name,
                                                               funcG2[region],funcN1[region],funcN2[region],
                                                               guessG2[region],guessN1[region],guessN2[region],
                                                               limits[region],region,show)
            # store values of fit parameters
            self.poptG2[region] = poptG2
            self.poptN1[region] = poptN1
            self.poptN2[region] = poptN2

        # add fit functions to class
        self.funcG2 = funcG2
        self.funcN1 = funcN1
        self.funcN2 = funcN2

    def ProdFitted(self,XY,region):
        """ Using the fit results of <N1><N2>, it computes <N1><N2>.
        Parameters
        --------------
        XY : numpy array
            axis values. XY should of the format [x,y] where x and y ate the x and y axis values.
        region : str
            region where <N1><N2> is to be computed. Options are "cross", "loc1" and "loc2"
        --------------
        Returns
        --------------
            <N1><N2> function
        """
        return self.funcProd[region](XY,*self.poptProd[region])
    
    def DCEN1(self,XY,region):
        """ Using the fit results, it computes <N1>.
        Parameters
        --------------
        XY : numpy array
            axis values. XY should of the format [x,y] where x and y ate the x and y axis values.
        region : str
            region where <N1> is to be computed. Options are "cross", "loc1" and "loc2"
        --------------
        Returns
        --------------
            <N1> function
        """
        return self.funcN1[region](XY,*self.poptN1[region])
    
    def DCEN2(self,XY,region):
        """ Using the fit results, it computes <N2>.
        Parameters
        --------------
        XY : numpy array
            axis values. XY should of the format [x,y] where x and y ate the x and y axis values.
        region : str
            region where <N2> is to be computed. Options are "cross", "loc1" and "loc2"
        --------------
        Returns
        --------------
            <N2> function
        """
        return self.funcN2[region](XY,*self.poptN2[region])
    
    def DCEProd(self,XY,region):
        """ Using the fit results for <N1> and <N2>, it computes <N1><N2>.
        Parameters
        --------------
        XY : numpy array
            axis values. XY should of the format [x,y] where x and y ate the x and y axis values.
        region : str
            region where <N1> and <N2> are to be computed. Options are "cross", "loc1" and "loc2"
        --------------
        Returns
        --------------
            <N1><N2> function
        """
        return self.DCEN1(XY,region)*self.DCEN2(XY,region)

    def G2Fitted(self,XY,region):
        """ Using the fit results, it computes G2.
        Parameters
        --------------
        XY : numpy array
            axis values. XY should of the format [x,y] where x and y ate the x and y axis values.
        region : str
            region where G2 is to be computed. Options are "cross", "loc1" and "loc2"
        --------------
        Returns
        --------------
            G2 function
        """
        return self.funcG2[region](XY,*self.poptG2[region])

    def g2Fitted(self,XY,region):
        """ Using the fit results of G2 and <N1><N2>, it computes g2.
        Parameters
        --------------
        XY : numpy array
            axis values. XY should of the format [x,y] where x and y ate the x and y axis values.
        region : str
            region where g2 is to be computed. Options are "cross", "loc1" and "loc2"
        --------------        
        Returns
        --------------
            g2 function
        """
        return self.G2Fitted(XY,region)/self.ProdFitted(XY,region)
    
    def g2DCE(self,XY,region):
        """ Using the fit results of G2, <N1> and <N2>, it computes g2.
        Parameters
        --------------
        XY : numpy array
            axis values. XY should of the format [x,y] where x and y ate the x and y axis values.
        region : str
            region where g2 is to be computed. Options are "cross", "loc1" and "loc2"
        --------------        
        Returns
        --------------
            g2 function
        """
        return self.G2Fitted(XY,region)/self.DCEProd(XY,region)

    def G2Integrated(self,axis,V0,U0,theta,width,region):
        """ Using the fit results, it integrates "G2" over the diagonal or anti-diagonal depending on the region chosen.
        If region is "cross" then we integrate over the diagonal, otherwise we integrate over the anti-diagonal.
        Parameters
        --------------
        axis : numpy array
            axis values: points along which integral is computed
        V0 , U0 : float
            center of integration region in rotated frame
        theta : float
            angle that diagonal does with x-axis
        width : float
            width of region of integration
        region : str
            region where integral is to be computed. Options are "cross", "loc1" and "loc2"
        --------------        
        Returns
        --------------
            G2 integrated integrated
        """

        # G2 integrated
        G2_1D = np.zeros(len(axis),dtype=float)

        # if region is cross, we integrate along the diagonal
        if region == "cross":
            # integrate along v = Vz1 + Vz2
            v = V0 + axis
            # for each value of axis
            for i in range(0,len(v)):
                # integrate G2
                IntFunc = lambda u : RotateFunction([v[i],u],theta,self.G2Fitted,[region])
                integral = quad(IntFunc , U0-width/2, U0+width/2)
                G2_1D[i] = integral[0]/width
        # else, we integrate along the anti-diagonal
        else:
            # integrate along u = Vz2 - Vz1
            u = U0 + axis
            # for each value of axis
            for i in range(0,len(u)):
                # integrate G2
                IntFunc = lambda v : RotateFunction([v,u[i]],theta,self.G2Fitted,[region])
                integral = quad(IntFunc , V0-width/2, V0+width/2)
                G2_1D[i] = integral[0]/width

        return G2_1D
    
    def ProdIntegrated(self,axis,V0,U0,theta,width,region):
        """ Using the fit results, it integrates "<N1><N2>" over the diagonal or anti-diagonal depending on the region chosen.
        If region is "cross" then we integrate over the diagonal, otherwise we integrate over the anti-diagonal.
        Parameters
        --------------
        axis : numpy array
            axis values: points along which integral is computed
        V0 , U0 : float
            center of integration region in rotated frame
        theta : float
            angle that diagonal does with x-axis
        width : float
            width of region of integration
        region : str
            region where integral is to be computed. Options are "cross", "loc1" and "loc2"
        --------------        
        Returns
        --------------
            <N1><N2> integrated integrated
        """

        # G2 integrated
        Prod_1D = np.zeros(len(axis),dtype=float)

        # if region is cross, we integrate along the diagonal
        if region == "cross":
            # integrate along v = Vz1 + Vz2
            v = V0 + axis
            # for each value of axis
            for i in range(0,len(v)):
                # integrate G2
                IntFunc = lambda u : RotateFunction([v[i],u],theta,self.ProdFitted,[region])
                integral = quad(IntFunc , U0-width/2, U0+width/2)
                Prod_1D[i] = integral[0]/width
        # else, we integrate along the anti-diagonal
        else:
            # integrate along u = Vz2 - Vz1
            u = U0 + axis
            # for each value of axis
            for i in range(0,len(u)):
                # integrate G2
                IntFunc = lambda v : RotateFunction([v,u[i]],theta,self.ProdFitted,[region])
                integral = quad(IntFunc , V0-width/2, V0+width/2)
                Prod_1D[i] = integral[0]/width

        return Prod_1D
    
    def DCEProdIntegrated(self,axis,V0,U0,theta,width,region):
        """ Using the fit results of <N1> and <N2>, it integrates "<N1><N2>" over the diagonal or anti-diagonal 
        depending on the region chosen.
        If region is "cross" then we integrate over the diagonal, otherwise we integrate over the anti-diagonal.
        Parameters
        --------------
        axis : numpy array
            axis values: points along which integral is computed
        V0 , U0 : float
            center of integration region in rotated frame
        theta : float
            angle that diagonal does with x-axis
        width : float
            width of region of integration
        region : str
            region where integral is to be computed. Options are "cross", "loc1" and "loc2"
        --------------        
        Returns
        --------------
            <N1><N2> integrated integrated
        """
        
        # <N1><N2> integrated
        Prod_1D = np.zeros(len(axis),dtype=float)

        # if region is cross, we integrate along the diagonal
        if region == "cross":
            # integrate along v = Vz1 + Vz2
            v = V0 + axis
            # for each value of axis
            for i in range(0,len(v)):
                # integrate G2
                IntFunc = lambda u : RotateFunction([v[i],u],theta,self.DCEProd,[region])
                integral = quad(IntFunc , U0-width/2, U0+width/2)
                Prod_1D[i] = integral[0]/width
        # else, we integrate along the anti-diagonal
        else:
            # integrate along u = Vz2 - Vz1
            u = U0 + axis
            # for each value of axis
            for i in range(0,len(u)):
                # integrate G2
                IntFunc = lambda v : RotateFunction([v,u[i]],theta,self.DCEProd,[region])
                integral = quad(IntFunc , V0-width/2, V0+width/2)
                Prod_1D[i] = integral[0]/width

        return Prod_1D
        
    def g2Integrated(self,axis,V0,U0,theta,width,region):
        """ Using the fit results, it integrates "g2" over the diagonal or anti-diagonal depending on the region chosen.
        If region is "cross" then we integrate over the diagonal, otherwise we integrate over the anti-diagonal.
        Parameters
        --------------
        axis : numpy array
            axis values: points along which integral is computed
        V0 , U0 : float
            center of integration region in rotated frame
        theta : float
            angle that diagonal does with x-axis
        width : float
            width of region of integration
        region : str
            region where integral is to be computed. Options are "cross", "loc1" and "loc2"
        --------------        
        Returns
        --------------
            g2 integrated integrated
        """

        # G2 integrated
        g2_1D = np.zeros(len(axis),dtype=float)

        # if region is cross, we integrate along the diagonal
        if region == "cross":
            # integrate along v = Vz1 + Vz2
            v = V0 + axis
            # for each value of axis
            for i in range(0,len(v)):
                # integrate G2
                IntFunc = lambda u : RotateFunction([v[i],u],theta,self.g2Fitted,[region])
                integral = quad(IntFunc , U0-width/2, U0+width/2)
                g2_1D[i] = integral[0]/width
        # else, we integrate along the anti-diagonal
        else:
            # integrate along u = Vz2 - Vz1
            u = U0 + axis
            # for each value of axis
            for i in range(0,len(u)):
                # integrate G2
                IntFunc = lambda v : RotateFunction([v,u[i]],theta,self.g2Fitted,[region])
                integral = quad(IntFunc , V0-width/2, V0+width/2)
                g2_1D[i] = integral[0]/width

        return g2_1D
    
    def g2DCEIntegrated(self,axis,V0,U0,theta,width,region):
        """ Using the fit results of <N1><N2>, it integrates "g2" over the diagonal or anti-diagonal depending on the region chosen.
        If region is "cross" then we integrate over the diagonal, otherwise we integrate over the anti-diagonal.
        Parameters
        --------------
        axis : numpy array
            axis values: points along which integral is computed
        V0 , U0 : float
            center of integration region in rotated frame
        theta : float
            angle that diagonal does with x-axis
        width : float
            width of region of integration
        region : str
            region where integral is to be computed. Options are "cross", "loc1" and "loc2"
        --------------        
        Returns
        --------------
            g2 integrated integrated
        """

        # G2 integrated
        g2_1D = np.zeros(len(axis),dtype=float)

        # if region is cross, we integrate along the diagonal
        if region == "cross":
            # integrate along v = Vz1 + Vz2
            v = V0 + axis
            # for each value of axis
            for i in range(0,len(v)):
                # integrate G2
                IntFunc = lambda u : RotateFunction([v[i],u],theta,self.g2DCE,[region])
                integral = quad(IntFunc , U0-width/2, U0+width/2)
                g2_1D[i] = integral[0]/width
        # else, we integrate along the anti-diagonal
        else:
            # integrate along u = Vz2 - Vz1
            u = U0 + axis
            # for each value of axis
            for i in range(0,len(u)):
                # integrate G2
                IntFunc = lambda v : RotateFunction([v,u[i]],theta,self.g2DCE,[region])
                integral = quad(IntFunc , V0-width/2, V0+width/2)
                g2_1D[i] = integral[0]/width

        return g2_1D
            
    def bootstrap_dataframes(self):
        """ bootstrap the dataframes [self.df_atoms_var1, self.df_atoms_var2]  in an efficient way.  """


        # if we do not have a copy of the 3 dataframes we create them before bootstrap
        if self.is_there_a_copy_of_total is False:
            # list to hold original cycles values
            self.cyclesList = []
            for j , (attr_atoms , attr_copy) in enumerate(zip(["df_atoms_var1", "df_atoms_var2"],["df_atoms_var1_copy", "df_atoms_var2_copy"])):
                # get atom dataframe ["cross", "loc1", "loc2"] from the object.
                df_atoms = getattr(self, attr_atoms)
                # create copy
                df_atoms_copy = copy.deepcopy(df_atoms)
                df_atoms_copy["Original Cycle"] = df_atoms_copy["Cycle"]
                # assign copy as atribute
                setattr(self, attr_copy,df_atoms_copy)
                self.cyclesList.append(copy.deepcopy(df_atoms["Cycle"]))
            self.is_there_a_copy_of_total = True
            print("[Warning] : I just saved a copy of the total dataframe because you will destruct your original dataframe.")


        """ bootstrap"""
        for j , (attr_atoms , attr_copy) in enumerate(zip(["df_atoms_var1", "df_atoms_var2"],["df_atoms_var1_copy", "df_atoms_var2_copy"])):   
            # get dataframe
            dfcopy =  getattr(self, attr_copy)  
            
            # copy dataset and groupby cycle
            ordata = dfcopy.set_index(["Cycle", dfcopy.groupby("Cycle").cumcount()])
            # get data in array format
            original_data_array = ordata.values.reshape((len(self.cyclesList[j]), -1, len(ordata.columns)))
            # get shape
            _, n_tmp_index, _ = original_data_array.shape
            # peak randomnly elements of the array
            new_indices = np.random.randint(0, len(self.cyclesList[j]), len(self.cyclesList[j]))
            NEW = np.zeros_like(original_data_array)
            NEW[:] = original_data_array[new_indices]
            # create new dataframe with randomized data
            df_atoms = pd.DataFrame(data=NEW.reshape((-1, len(ordata.columns))), columns=list(ordata.columns))
            df_atoms["Cycle"] = np.repeat(self.cyclesList[j], n_tmp_index)

            # assign dataframe as atribute
            setattr(self, attr_atoms,df_atoms)  

    def bootstrap_dataframe_atoms(self):
        """ bootstrap the dataframes self.atoms  in an efficient way.  """

        # if we do not have a copy of the 3 dataframes we create them before bootstrap
        if self.is_there_a_copy_of_total is False:
            # create copy
            self.atoms_copy = copy.deepcopy(self.atoms)
            self.atoms_copy["Original Cycle"] = self.atoms_copy["Cycle"]
            self.cyclesList = copy.deepcopy(self.atoms_copy["Cycle"].drop_duplicates().to_list())
            self.is_there_a_copy_of_total = True
            print("[Warning] : I just saved a copy of the total dataframe because you will destruct your original dataframe.")


        """ bootstrap"""    
        # randomnly select cycles
        newCycles = np.random.choice(self.cyclesList, size = len(self.cyclesList))    
        # create temporary dataset to preserve the index by assigning it as a column first so we can set_index after the merging
        dfTemp = pd.DataFrame({"Cycle": newCycles, "OCycles": self.cyclesList})
        # merge dataset "dfTemp" with "self.atoms_copy" dataframe so that we keep only cycles in "newCycles" and preserve the size of the dataframe
        self.atoms =  pd.merge(dfTemp, self.atoms_copy, on="Cycle").reset_index().drop(["index","Original Cycle","Cycle"], axis=1).rename(columns={'OCycles':'Cycle'})


# returns a rotated function
def RotateFunction(VU,theta,func,parmaters):
    """ computes the the function func, defined in cartesisan axis x and y, in the rotated axis
    v and u. """
    v = VU[0]
    u = VU[1]
    # rotate axis
    x , y = rotation(v,u,-theta)
    return func([x,y],*parmaters)

# return a rotated g2
def g2Rotated(VU,theta,poptProd,poptG2,funcG2,funcProd):
    """ given two fit functions, poptProd and funG2, and their fit parameters
    it rotates the axis and returns another a rotated g2"""
    v = VU[0]
    u = VU[1]
    # compute <N1><N2>
    Prod = RotateFunction([v,u],theta,funcProd,poptProd)
    # compute G2
    G2 = RotateFunction([v,u],theta,funcG2,poptG2)
    # compute normalized g2    
    g2 = G2/Prod
    return g2


# function to fit G2 and <N1> and <N2> at the same time
def FitG2andDCEPeaks(df,xname,yname,funcG2,funcN1,funcN2,guessG2,guessN1,guessN2,limits,title = "title",show = False):
    """ Fits a function funcG2 to functions G2, <N1> and <N2> at the same time to ensure that G2 and
    <N1><N2> have the same offset
    Parameters
    --------------
    df : pandas dataframe
        dataframe with data to fit
    xname : str
        name of collumn with x-axis values
    yname : str
        name of collumn with y-axis values
    funcG2 : function used to fit G2
        funcG2 should be a function of two variables, x and y, that should be passed as an array [x,y]. 
        Rest of the arguments for funcG2 are considered to be fit parameters. E.g: funcG2 can be a Skew2DGaussian
    funcN1 : function used to fit <N1>
        Same type of function as funcG2. E.g: func can be a DCE_N1
    funcN2 : function used to fit <N2>
        Same type of function as funcG2. E.g: func can be a DCE_N2
    guessG2 : 1D numpy array
        initial guess for fit parameters of G2, offset included
    guessN1 : 1D numpy array
        initial guess for fit parameters of <N1>, offset included
    guessN2 : 1D numpy array
        initial guess for fit parameters of <N2>, offset included
    limits : 2D numpy array
        limits is of the form [limit_x,limit_y]. limit_x are the limits of the region of interest of the fit on the x-axis. 
        Equivalent for limit_y on the y-axis.
    region : str
        region where fit is to be performed: "cross", "loc1" or "loc2"
    title : str
        title of plot
    show : bool
        True if you want to plot fit results. False otherwise.
    Returns
    --------------
    poptG2 : 1D numpy array
        optimized fit parameters found by the fit for G2
    poptN1 : 1D numpy array
        optimized fit parameters found by the fit for <N1>
    poptN2 : 1D numpy array
        optimized fit parameters found by the fit for <N2>
    pcov : 2D numpy array
        covariance matrix of fit parameters (both G2 and <N1> and <N2> )
    """

    # name of quantities to fit
    keys = ["G2","N_1","N_2"]

    # initialize arrays to fit
    comboAxis = [[],[]]
    comboData = []

    # apply limits to data
    limit_x = limits[0]
    limit_y = limits[1]
    dfFilter = df.loc[(df[xname] > limit_x[0]) & (df[xname] < limit_x[1]) & (df[yname] > limit_y[0]) & (df[yname] < limit_y[1]) ]
        
    # for eack key get data
    for j in range(0,len(keys)):
        # get axis and data values for G2
        xdata , ydata , zdata , xy_stack , zflat = FlattendData(dfFilter,keys[j],xname,yname)
        # concatenate data to fit
        comboAxis = np.concatenate((comboAxis,xy_stack),axis=1)
        comboData = np.concatenate((comboData,zflat))

    
    # number of G2 fit parameters
    NG2 = len(guessG2) 
    # number of <N1> fit parameters
    NN1 = len(guessN1)
    # number of <N2> fit parameters
    NN2 = len(guessN2)
    # create Fit function object 
    FitObject = fitClass_v2(funcG2,funcN1,funcN2,NG2,NN1,NN2)

    # build guess array for combined fit. We ignore the offset in guessG2
    guess = np.concatenate((guessG2[0:NG2-1],guessN1,guessN2))

    # try to fit data
    try:
        popt , pcov = curve_fit(FitObject.combinedFunction,comboAxis,comboData,p0=guess)
    except:
        print("Fit failed: Couldn't fit both at the same time!")
        popt = guess
        pcov = []

    # seperate parameters of G2, <N1> and <N2>
    poptG2 = popt[:NG2-1]
    poptN1 = popt[NG2-1:NG2-1+NN1]
    poptN2 = popt[NG2-1+NN1:]

    # add offset to G2 parameters
    offset = poptN1[-1]*poptN2[-1]
    poptG2 = np.concatenate((poptG2,[offset]))
          
    # if user wants plot
    if show:
        fig, ax = plt.subplots(nrows = len(keys),ncols = 4,figsize = (20,20))

        # for eack key get data
        func = [funcG2,funcN1,funcN2]
        popt = [poptG2,poptN1,poptN2]

        """ plot G2 """
        if True:
            # get axis and data values for G2
            xdata , ydata , zdata , xy_stack , zflat = FlattendData(df,keys[0],xname,yname)
            """ plot original data """
            X , Y = np.meshgrid(xdata, ydata)
            plot1 = ax[0,0].pcolormesh(X, Y, zdata)
            ax[0,0].set_xlabel('Vz1')
            ax[0,0].set_ylabel('Vz2')
            ax[0,0].set_title(keys[0]+": data")
            cb1 = fig.colorbar(plot1,ax=ax[0,0])
            ax[0,0].contour(X, Y, zdata,levels = 3,colors = "white",linestyles = "dashed",linewidths = 0.8)
            ax[0,0].set_title(keys[0])
            """ plot fit result """
            # compute fit function
            x2 = np.linspace(np.amin(xdata),np.amax(xdata),500)
            y2 = np.linspace(np.amin(ydata),np.amax(ydata),500)
            X1, X2 = np.meshgrid(x2, y2)
            Z = func[0]([X1,X2],*popt[0])
            # plot
            plot2 = ax[0,1].pcolormesh(X1, X2, Z)
            ax[0,1].set_xlabel('Vz1')
            ax[0,1].set_title(keys[0]+": fit")
            cb2 = fig.colorbar(plot2,ax=ax[0,1])
            cb2.mappable.set_clim(*cb1.mappable.get_clim())
            ax[0,1].contour(X1, X2, Z,levels = 3,colors = "white",linestyles = "dashed",linewidths = 0.8)
            """ plot interpolation of data """
            # compute 2D interpolation
            funcInter = RegularGridInterpolator((xdata,ydata),np.transpose(zdata),method = "cubic")
            Z = funcInter((X1,X2))
            # plot
            plot3 = ax[0,2].pcolormesh(X1, X2, Z)
            ax[0,2].set_xlabel('Vz1')
            ax[0,2].set_title(keys[0]+": interpolation")
            cb3 = fig.colorbar(plot3,ax=ax[0,2])
            cb3.mappable.set_clim(*cb1.mappable.get_clim())
            ax[0,2].contour(X1, X2, Z,levels = 3,colors = "white",linestyles = "dashed",linewidths = 0.8)
            """ plot 2 sigma region to help plot visualization """
            amin = poptN1[0] - 2*poptN1[1]
            amax = poptN1[0] + 2*poptN1[1]
            bmin = poptN2[0] - 2*poptN2[1]
            bmax = poptN2[0] + 2*poptN2[1]
            for i in [0,1,2]:
                ax[0,i].hlines(bmin,xmin=amin,xmax=amax,color = "red",ls = "dashed")
                ax[0,i].hlines(bmax,xmin=amin,xmax=amax,color = "red",ls = "dashed")
                ax[0,i].vlines(amin,ymin=bmin,ymax=bmax,color = "red",ls = "dashed")
                ax[0,i].vlines(amax,ymin=bmin,ymax=bmax,color = "red",ls = "dashed")
            """ plot difference between data and fit """
            Deltaz = np.abs(func[0]([X,Y],*popt[0]) - zdata)/func[0]([X,Y],*popt[0])
            plot4 = ax[0,3].pcolormesh(X, Y, Deltaz,vmin = 0.0,vmax = 0.25)
            ax[0,3].set_xlabel('Vz1')
            ax[0,3].set_title(keys[0]+": (data - fit)/fit")
            cb4 = fig.colorbar(plot4,ax=ax[0,3])
            ax[0,3].contour(X, Y, zdata,levels = 3,colors = "white",linestyles = "dashed",linewidths = 0.8)

            # fix axis limits
            for i in [0,1,2,3]:
                ax[0,i].set_xlim(np.amin(xdata),np.amax(xdata))
                ax[0,i].set_ylim(np.amin(ydata),np.amax(ydata))
            
        """ plot <N1>  """
        if True:
            # get axis and data values for <N1>
            xdata , ydata , zdata , xy_stack , zflat = FlattendData(df,keys[1],xname,yname)
            """ plot original data """
            ax[1,0].plot(xdata, zdata[0],label = "data")
            """ plot fit result """
            # compute fit function
            x2 = np.linspace(np.amin(xdata),np.amax(xdata),500)
            y2 = np.linspace(np.amin(ydata),np.amax(ydata),500)
            Z = func[1]([x2,y2[0]],*popt[1])
            # plot
            ax[1,0].plot(x2, Z,label = "fit")
            ax[1,0].set_xlabel('Vz')
            ax[1,0].set_ylabel('Density')
            ax[1,0].set_title(keys[1])
            ax[1,0].legend()

            """ plot 2 sigma region to help plot visualization """
            amin = poptN1[0] - 2*poptN1[1]
            amax = poptN1[0] + 2*poptN1[1]
            ax[0,0].vlines(amin,ymin=0,ymax=np.amax(Z),color = "red",ls = "dashed")
            ax[0,0].vlines(amax,ymin=0,ymax=np.amax(Z),color = "red",ls = "dashed")

            """ plot difference between data and fit """
            Deltaz = np.abs(func[1]([xdata,ydata[0]],*popt[1]) - zdata[0])/func[1]([xdata,ydata[0]],*popt[1])
            ax[1,1].plot(xdata,Deltaz)
            ax[1,1].set_xlabel('Vz')
            ax[1,1].set_title(keys[1]+": (data - fit)/fit")
            ax[1,1].vlines(amin,ymin=0,ymax=1,color = "red",ls = "dashed")
            ax[1,1].vlines(amax,ymin=0,ymax=1,color = "red",ls = "dashed")
            ax[1,1].set_ylim(0,1)

        """ plot <N2> """
        if True:
            # get axis and data values for <N2>
            xdata , ydata , zdata , xy_stack , zflat = FlattendData(df,keys[2],xname,yname)
            """ plot original data """
            ax[1,2].plot(ydata, zdata[:,0],label = "data")
            """ plot fit result """
            # compute fit function
            x2 = np.linspace(np.amin(xdata),np.amax(xdata),500)
            y2 = np.linspace(np.amin(ydata),np.amax(ydata),500)
            Z = func[2]([x2[0],y2],*popt[2])
            # plot
            ax[1,2].plot(y2, Z,label = "fit")
            ax[1,2].legend()
            ax[1,2].set_xlabel('Vz')
            ax[1,2].set_ylabel('Density')
            ax[1,2].set_title(keys[2])

            """ plot 2 sigma region to help plot visualization """
            amin = poptN2[0] - 2*poptN2[1]
            amax = poptN2[0] + 2*poptN2[1]
            ax[1,2].vlines(amin,ymin=0,ymax=np.amax(Z),color = "red",ls = "dashed")
            ax[1,2].vlines(amax,ymin=0,ymax=np.amax(Z),color = "red",ls = "dashed")

            """ plot difference between data and fit """
            Deltaz = np.abs(func[2]([xdata[0],ydata],*popt[2]) - zdata[:,0])/func[2]([xdata[0],ydata],*popt[2])
            ax[1,3].plot(ydata,Deltaz)
            ax[1,3].set_xlabel('Vz')
            ax[1,3].set_title(keys[2]+": (data - fit)/fit")
            ax[1,3].vlines(amin,ymin=0,ymax=1,color = "red",ls = "dashed")
            ax[1,3].vlines(amax,ymin=0,ymax=1,color = "red",ls = "dashed")
            ax[1,3].set_ylim(0,1)

        """ plot product <N1><N2> as well"""
        if True:
            # get axis and data values for G2
            xdata , ydata , zdata , xy_stack , zflat = FlattendData(df,"<N1><N2>",xname,yname)
            j = 2
            """ plot original data """
            X , Y = np.meshgrid(xdata, ydata)
            plot1 = ax[j,0].pcolormesh(X, Y, zdata)
            ax[j,0].set_xlabel('Vz1')
            ax[j,0].set_ylabel('Vz2')
            ax[j,0].set_title("<N1><N2>: data")
            cb1 = fig.colorbar(plot1,ax=ax[j,0])
            ax[j,0].contour(X, Y, zdata,levels = 3,colors = "white",linestyles = "dashed",linewidths = 0.8)
            ax[j,0].set_title("<N1><N2>")

            """ plot fit result """
            # compute fit function
            x2 = np.linspace(np.amin(xdata),np.amax(xdata),500)
            y2 = np.linspace(np.amin(ydata),np.amax(ydata),500)
            X1, X2 = np.meshgrid(x2, y2)
            Z = funcN1([X1,X2],*poptN1)*funcN2([X1,X2],*poptN2)
            # plot
            plot2 = ax[j,1].pcolormesh(X1, X2, Z)
            ax[j,1].set_xlabel('Vz1')
            ax[j,1].set_title("<N1><N2>")
            cb2 = fig.colorbar(plot2,ax=ax[j,1])
            cb2.mappable.set_clim(*cb1.mappable.get_clim())
            ax[j,1].contour(X1, X2, Z,levels = 3,colors = "white",linestyles = "dashed",linewidths = 0.8)

            """ plot interpolation of data """
            # compute 2D interpolation
            funcInter = RegularGridInterpolator((xdata,ydata),np.transpose(zdata),method = "cubic")
            Z = funcInter((X1,X2))
            # plot
            plot3 = ax[j,2].pcolormesh(X1, X2, Z)
            ax[j,2].set_xlabel('Vz1')
            ax[j,2].set_title("<N1><N2>: interpolation")
            cb3 = fig.colorbar(plot3,ax=ax[j,2])
            cb3.mappable.set_clim(*cb1.mappable.get_clim())
            ax[j,2].contour(X1, X2, Z,levels = 3,colors = "white",linestyles = "dashed",linewidths = 0.8)

            """ plot 2 sigma region to help plot visualization """
            amin = poptN1[0] - 2*poptN1[1]
            amax = poptN1[0] + 2*poptN1[1]
            bmin = poptN2[0] - 2*poptN2[1]
            bmax = poptN2[0] + 2*poptN2[1]
            for i in [0,1,2]:
                ax[j,i].hlines(bmin,xmin=amin,xmax=amax,color = "red",ls = "dashed")
                ax[j,i].hlines(bmax,xmin=amin,xmax=amax,color = "red",ls = "dashed")
                ax[j,i].vlines(amin,ymin=bmin,ymax=bmax,color = "red",ls = "dashed")
                ax[j,i].vlines(amax,ymin=bmin,ymax=bmax,color = "red",ls = "dashed")

            """ plot difference between data and fit """
            Zfit = funcN1([X,Y],*poptN1)*funcN2([X,Y],*poptN2)
            Deltaz = np.abs(Zfit - zdata)/Zfit
            plot4 = ax[j,3].pcolormesh(X, Y, Deltaz,vmin = 0.0,vmax = 0.25)
            ax[j,3].set_xlabel('Vz1')
            ax[j,3].set_title("<N1><N2>: (data - fit)/fit")
            cb4 = fig.colorbar(plot4,ax=ax[j,3])
            ax[j,3].contour(X, Y, zdata,levels = 3,colors = "white",linestyles = "dashed",linewidths = 0.8)
        
        fig.suptitle(title)
        plt.show()

    return poptG2 , poptN1 , poptN2 , pcov

# function to fit G2 and <N1><N2> at the same time
def FitG2andProd(df,xname,yname,funcG2,funcProd,guessG2,guessProd,limits,title = "title",show = False):
    """ Fits a function funcG2 to functions G2 and funcProd <N1><N2> at the same time to ensure that they have the same offset
    Parameters
    --------------
    df : pandas dataframe
        dataframe with data to fit
    xname : str
        name of collumn with x-axis values
    yname : str
        name of collumn with y-axis values
    funcG2 : function used to fit G2
        funcG2 should be a function of two variables, x and y, that should be passed as an array [x,y]. 
        Rest of the arguments for funcG2 are considered to be fit parameters. E.g: funcG2 can be a Skew2DGaussian
    funcProd : function used to fit <N1><N2>
        Same type of function as funcG2. E.g: func can be a ProdGaussian
    guessG2 : 1D numpy array
        initial guess for fit parameters of G2, offset included
    guessProd : 1D numpy array
        initial guess for fit parameters of <N1><N2>, offset included
    limits : 2D numpy array
        limits is of the form [limit_x,limit_y]. limit_x are the limits of the region of interest of the fit on the x-axis. 
        Equivalent for limit_y on the y-axis.
    title : str
        title of plot
    show : bool
        True if you want to plot fit results. False otherwise.
    Returns
    --------------
    poptG2 : 1D numpy array
        optimized fit parameters found by the fit for G2
    poptProd : 1D numpy array
        optimized fit parameters found by the fit for <N1><N2>
    pcov : 2D numpy array
        covariance matrix of fit parameters (both G2 and <N1><N2> )
    """

    # name of quantities to fit
    keys = ["G2","<N1><N2>"]

    # initialize arrays to fit
    comboAxis = [[],[]]
    comboData = []

    # apply limits to data
    limit_x = limits[0]
    limit_y = limits[1]
    dfFilter = df.loc[(df[xname] > limit_x[0]) & (df[xname] < limit_x[1]) & (df[yname] > limit_y[0]) & (df[yname] < limit_y[1]) ]
        
    # for eack key get data
    for j in range(0,len(keys)):
        # get axis and data values for G2
        xdata , ydata , zdata , xy_stack , zflat = FlattendData(dfFilter,keys[j],xname,yname)
        # concatenate data to fit
        comboAxis = np.concatenate((comboAxis,xy_stack),axis=1)
        comboData = np.concatenate((comboData,zflat))

    
    # number of G2 fit parameters
    NG2 = len(guessG2) 
    # number of <N1><N2> fit parameters
    NProd = len(guessProd)
    # create Fit function object 
    FitObject = fitClass(funcG2,funcProd,NG2,NProd)    
    # build guess array for combined fit
    guess = np.concatenate((guessG2[0:NG2-1],guessProd))
    # try to fit data
    try:
        popt , pcov = curve_fit(FitObject.combinedFunction,comboAxis,comboData,p0=guess)
    except:
        print("Fit failed: Couldn't fit both at the same time!")
        popt = guess
        pcov = []

    # seperate parameters of G2 from <N1><N2>
    poptG2 = np.concatenate((popt[:int(NG2-1)],[popt[-1]]))
    poptProd = popt[int(NG2-1):]
          
    # if user wants plot
    if show:
        fig, ax = plt.subplots(nrows=len(keys),ncols = 4,sharey = True,sharex=True,figsize = (20,9))

        # for eack key get data
        func = [funcG2,funcProd]
        popt = [poptG2,poptProd]
        for j in range(0,len(keys)):
            # get axis and data values for G2
            xdata , ydata , zdata , xy_stack , zflat = FlattendData(df,keys[j],xname,yname)

            """ plot original data """
            X , Y = np.meshgrid(xdata, ydata)
            plot1 = ax[j,0].pcolormesh(X, Y, zdata)
            ax[j,0].set_xlabel('Vz1')
            ax[j,0].set_ylabel('Vz2')
            ax[j,0].set_title(keys[j]+": data")
            cb1 = fig.colorbar(plot1,ax=ax[j,0])
            ax[j,0].contour(X, Y, zdata,levels = 3,colors = "white",linestyles = "dashed",linewidths = 0.8)
            ax[j,0].set_title(keys[j])

            """ plot fit result """
            # compute fit function
            x2 = np.linspace(np.amin(xdata),np.amax(xdata),500)
            y2 = np.linspace(np.amin(ydata),np.amax(ydata),500)
            X1, X2 = np.meshgrid(x2, y2)
            Z = func[j]([X1,X2],*popt[j])
            # plot
            plot2 = ax[j,1].pcolormesh(X1, X2, Z)
            ax[j,1].set_xlabel('Vz1')
            ax[j,1].set_title(keys[j]+": fit")
            cb2 = fig.colorbar(plot2,ax=ax[j,1])
            cb2.mappable.set_clim(*cb1.mappable.get_clim())
            ax[j,1].contour(X1, X2, Z,levels = 3,colors = "white",linestyles = "dashed",linewidths = 0.8)

            """ plot interpolation of data """
            # compute 2D interpolation
            funcInter = RegularGridInterpolator((xdata,ydata),np.transpose(zdata),method = "cubic")
            Z = funcInter((X1,X2))
            # plot
            plot3 = ax[j,2].pcolormesh(X1, X2, Z)
            ax[j,2].set_xlabel('Vz1')
            ax[j,2].set_title(keys[j]+": interpolation")
            cb3 = fig.colorbar(plot3,ax=ax[j,2])
            cb3.mappable.set_clim(*cb1.mappable.get_clim())
            ax[j,2].contour(X1, X2, Z,levels = 3,colors = "white",linestyles = "dashed",linewidths = 0.8)

            """ plot 2 sigma region to help plot visualization """
            amin = poptProd[0] - 2*poptProd[3]
            amax = poptProd[0] + 2*poptProd[3]
            bmin = poptProd[1] - 2*poptProd[4]
            bmax = poptProd[1] + 2*poptProd[4]
            for i in [0,1,2]:
                ax[j,i].hlines(bmin,xmin=amin,xmax=amax,color = "red",ls = "dashed")
                ax[j,i].hlines(bmax,xmin=amin,xmax=amax,color = "red",ls = "dashed")
                ax[j,i].vlines(amin,ymin=bmin,ymax=bmax,color = "red",ls = "dashed")
                ax[j,i].vlines(amax,ymin=bmin,ymax=bmax,color = "red",ls = "dashed")

            """ plot difference between data and fit """
            Deltaz = np.abs(func[j]([X,Y],*popt[j]) - zdata)/func[j]([X,Y],*popt[j])
            plot4 = ax[j,3].pcolormesh(X, Y, Deltaz,vmin = 0.0,vmax = 0.25)
            ax[j,3].set_xlabel('Vz1')
            ax[j,3].set_title(keys[j]+": (data - fit)/fit")
            cb4 = fig.colorbar(plot4,ax=ax[j,3])
            ax[j,3].contour(X, Y, zdata,levels = 3,colors = "white",linestyles = "dashed",linewidths = 0.8)
        
        fig.suptitle(title)
        plt.show()

    return poptG2 , poptProd , pcov

# Another verison of fit class: defined so that the code can be used in a more dynamic way
class fitClass_v2:
    """
    Parameters
    ----------
    funcG2 : function used to fit G2
        func should be a function of two variables, x and y, that should be passed as an array [x,y]. 
        Rest of the arguments for func are considered to be fit parameters and offset parameter should be the last parameter
        e.g: funcG2 can be a Skew2DGaussian
    funcN1 : function used to fit <N1>
        same type of function as funcG2. E.g: funcProd can be a DCE_N1
    funcN2: function used to fit <N2>
        same type of function as funcG2. E.g: funcProd can be a DCE_N2
    NG2 : int
        number of fit parameters to be passed to function funcG2
    NN1 : int
        number of fit parameters to be passed to function funcN1
    NN2 : int
        number of fit parameters to be passed to function funcN2
    """

    # initialize class object
    def __init__(self,funcG2,funcN1,funcN2,NG2,NN1,NN2):
        self.funcG2 = funcG2
        self.funcN1 = funcN1
        self.funcN2 = funcN2
        self.NG2 = int(NG2)
        self.NN1 = int(NN1)
        self.NN2 = int(NN2)

    # function used to fit G2, <N1> and <N2> at the same time so G2 and <N1><N2> have the same offset
    def combinedFunction(self,comboAxis,*fitParameters):
        """ function used to fit both G2, <N1> and <N2> at the same time so G2 and <N1><N2> have the same offset.
        G2, <N1> and <N2> can be fitted using different function, 
        but the user should make sure that all functions have offset parameter as last one.
        Parameters
        --------------
        comboAxis : 1D numpy array
            array with concatenated of two individual arrays with values of axis for G2, <N1> and <N2>. 
            G2 values should be first, followed by <N1> and finally <N2>.
            Each individual array should be a stack of x and y-axis values and it is assumed that they have equal length.
        fitParameters : list
            list of fit parameters passed on to funcG2 to fit G2 data and on to funcN1 and funcN2 to fit <N1> and <N2> data.
            G2 parameters, except the offset, should be first, followed by the <N1> paramaters and then finnaly by <N2> parameters.
            G2 offset will be the product of <N1> and <N2> offsets.  
        Returns
        --------------
            Concatenated array with G2, <N1> and <N2> values. G2 values are first.
        """

        # pack input parameters
        par = [*fitParameters]

        # get fit parameters for G2
        fitG2 = par[:self.NG2-1]

        # get fit parameters for <N1>
        fitN1 = par[self.NG2-1:self.NG2-1+self.NN1]
        # get fit parameters for <N2>
        fitN2 = par[self.NG2-1+self.NN1:]

        # add offset to G2 parameters
        offset = fitN1[-1]*fitN2[-1]
        fitG2 = np.concatenate((fitG2,[offset])) # add offset

        # get length of individual arrays
        length = int(np.shape(comboAxis)[1]/3)
        # get stacked axis-values for G2
        xy_G2 = comboAxis[:,:length]
        # get stacked axis-values for <N1> 
        xy_N1 = comboAxis[:,length:2*length]
        # get stacked axis-values for <N2> 
        xy_N2 = comboAxis[:,2*length:]
        
        # calculate flatten G2 array 
        G2 = self.funcG2(xy_G2,*fitG2)

        # calculate flatten <N1> array
        N1 = self.funcN1(xy_N1,*fitN1)
        # calculate flatten <N2> array
        N2 = self.funcN2(xy_N2,*fitN2)

        # return concatenated function values
        return np.concatenate((G2,N1,N2))

# Fit class: defined so that the code can be used in a more dynamic way
class fitClass:
    """
    Parameters
    ----------
    funcG2 : function used to fit G2
        func should be a function of two variables, x and y, that should be passed as an array [x,y]. 
        Rest of the arguments for func are considered to be fit parameters and offset parameter should be the last parameter
        e.g: funcG2 can be a Skew2DGaussian
    funcProd : function used to fit <N1><N2>
        same type of function as funcG2. E.g: funcProd can be a ProdGaussian
    NG2 : int
        number of fit parameters to be passed to function funcG2
    NProd : int
        number of fit parameters to be passed to function funcProd
    """

    # initialize class object
    def __init__(self,funcG2,funcProd,NG2,NProd):
        self.funcG2 = funcG2
        self.funcProd = funcProd
        self.NG2 = NG2
        self.Nprod = NProd

    # function used to fit G2 and <N1><N2> at the same time so they have the same offset
    def combinedFunction(self,comboAxis,*fitParameters):
        """ function used to fit both G2 and <N1><N2> at the same time so they have the same offset. G2 and <N1><N2> can
        be fitted using different function, but the user should make sure that both functions have equal offset limit.
        Parameters
        --------------
        comboAxis : 1D numpy array
            array with concatenated of two individual arrays with values of axis for G2 and <N1><N2>. G2 values should be first.
            Each individual array should be a stack of x and y-axis values and it is assumed that they have equal length.
        fitParameters : list
            list of fit parameters passed on to funcG2 to fit G2 data and on to funcProd to fit <N1><N2> data.
            G2 parameters, except the offset, should be first, followed by the <N1><N2> paramaters and finally the offset which will be the same for both.
        Returns
        --------------
            Concatenated array with G2 and <N1><N2> values. G2 values are first.
        """

        # pack input parameters
        par = [*fitParameters]

        # get fit parameters for G2
        fitG2 = par[:int(self.NG2-1)]
        fitG2 = np.concatenate((fitG2,[par[-1]])) # add offset

        # get fit parameters for <N1><N2>
        fitProd = par[int(self.NG2-1):]

        # get length of individual arrays
        length = int(np.shape(comboAxis)[1]/2)
        # get stacked axis-values for G2
        xy_G2 = comboAxis[:,:length]
        # get stacked axis-values for <N1><N2>
        xy_Prod = comboAxis[:,length:]

        # calculate flatten G2 array 
        G2 = self.funcG2(xy_G2,*fitG2)

        # calculate flatten <N1><N2> array
        Prod = self.funcProd(xy_Prod,*fitProd)

        # return concatenated function values
        return np.concatenate((G2,Prod))
    
# function that flattens data so that we can fit a 2D function to it
def FlattendData(df,key,xname,yname):
    """ Transforms dataframe collumns into numpy arrays and fllatens them so they can be used to fit a 2D fucntion.
    Parameters
    --------------
    df : pandas dataframe
        dataframe with data to fit
    key : str
        name of the collumn with data to fit. This is the z-axis of the function
    xname : str
        name of collumn with x-axis values
    yname : str
        name of collumn with y-axis values
    Returns
    --------------
    xdata : 1D numpy array
        x-axis values
    ydata : 1D numpy array
        y-axis values
    zdata : 2D numpy array
        data values
    xy_stack: stacked numpy array
        flattened and stacked array with x-axis and y-axis values for fit
    z_1D : 1D numpy array
        flattened array with z-axis values for fit
    """

    # from dataframe get 2D table of key data with rows = yname and collumns = xname
    quant = df.pivot(index=yname, columns=xname, values=key)
    # transform table to 2D numpy array which is now our 2D data
    zdata = quant.to_numpy()
    # get x-axis values
    ydata = np.array(quant.index)
    # get y-axis values
    xdata = np.array(quant.keys())

    # create mesh grid from x  and y axis values
    X , Y = np.meshgrid(xdata, ydata)
    # get shape of grid
    size = X.shape
    # resize X and Y to flatten them
    x_1d = X.reshape((1, np.prod(size)))
    y_1d = Y.reshape((1,np.prod(size)))
    # stack axis
    xy_stack = np.vstack((x_1d, y_1d))
    # flatten data to fit
    z_1D = zdata.flatten()

    return xdata , ydata , zdata , xy_stack , z_1D

# function to fit 2D data
def fit2D(func,df,key,xname,yname,guess,limits,title = "title",show = False):
    """ Given a scalar function func of 2 varibales, fits it to 2D data
    Parameters
    --------------
    func : function passed to curve_fit
        func should be a function of two variables, x and y, that should be passed as an array [x,y]. Rest of the arguments for func are considered to be fit parameters
    df : pandas dataframe
        dataframe with data to fit
    key : str
        name of the collumn with data to fit
    xname : str
        name of collumn with x-axis values
    yname : str
        name of collumn with y-axis values
    guess : 1D numpy array
        initial guess for fit parameters
    limits : 2D numpy array
        limits is of the form [limit_x,limit_y]. limit_x are the limits of the region of interest of the fit on the x-axis. Equivalent for limit_y on the y-axis.
    show : bool
        True if you want to plot fit results. False otherwise.
    Returns
    --------------
    popt : 1D numpy array
        optimized fit parameters found by the fit
    pcov : 2D numpy array
        covariance matrix of fit parameters
    """

    # apply limits to data
    limit_x = limits[0]
    limit_y = limits[1]
    dfFilter = df.loc[(df[xname] > limit_x[0]) & (df[xname] < limit_x[1]) & (df[yname] > limit_y[0]) & (df[yname] < limit_y[1]) ]

    # get axis and data values and flatten them to fit 2D function
    xdata , ydata , zdata , xy_stack , zflat = FlattendData(dfFilter,key,xname,yname)
    
    # try to fit
    try:
        popt , pcov = curve_fit(func,xy_stack,zflat,p0=guess)
    except:
        print("Fit failed !!")
        popt = guess
        pcov = []
    
    # plot data, fit, interpolation and data-fit
    if show:
        # get axis and data values and flatten them to fit 2D function
        xdata , ydata , zdata , xy_stack , zflat = FlattendData(df,key,xname,yname)

        fig, ax = plt.subplots(ncols = 4,sharey = True,sharex=True,figsize = (20,4))

        """ plot original data """
        X , Y = np.meshgrid(xdata, ydata)
        plot1 = ax[0].pcolormesh(X, Y, zdata)
        ax[0].set_xlabel('Vz1')
        ax[0].set_ylabel('Vz2')
        ax[0].set_title("data")
        cb1 = fig.colorbar(plot1,ax=ax[0])
        ax[0].contour(X, Y, zdata,levels = 3,colors = "white",linestyles = "dashed",linewidths = 0.8)

        """ plot fit result """

        # compute fit function
        x2 = np.linspace(np.amin(xdata),np.amax(xdata),500)
        y2 = np.linspace(np.amin(ydata),np.amax(ydata),500)
        X1, X2 = np.meshgrid(x2, y2)
        Z = func([X1,X2],*popt)

        plot2 = ax[1].pcolormesh(X1, X2, Z)
        ax[1].set_xlabel('Vz1')
        ax[1].set_title("fit")
        cb2 = fig.colorbar(plot2,ax=ax[1])
        cb2.mappable.set_clim(*cb1.mappable.get_clim())
        ax[1].contour(X1, X2, Z,levels = 3,colors = "white",linestyles = "dashed",linewidths = 0.8)

        """ plot interpolation of data """

        # compute 2D interpolation
        funcInter = RegularGridInterpolator((xdata,ydata),np.transpose(zdata),method = "cubic")
        Z = funcInter((X1,X2))

        plot3 = ax[2].pcolormesh(X1, X2, Z)
        ax[2].set_xlabel('Vz1')
        ax[2].set_title("interpolation")
        cb3 = fig.colorbar(plot3,ax=ax[2])
        cb3.mappable.set_clim(*cb1.mappable.get_clim())
        ax[2].contour(X1, X2, Z,levels = 3,colors = "white",linestyles = "dashed",linewidths = 0.8)

        """ plot difference between data and fit """
        Deltaz = np.abs(func([X,Y],*popt) - zdata)/func([X,Y],*popt)
        plot4 = ax[3].pcolormesh(X, Y, Deltaz,vmin = 0.0,vmax = 0.25)
        ax[3].set_xlabel('Vz1')
        ax[3].set_title("(data - fit)/fit")
        #ax[3].set_xlim(-14,-10)
        #ax[3].set_ylim(10,14)
        cb4 = fig.colorbar(plot4,ax=ax[3])
        ax[3].contour(X, Y, zdata,levels = 3,colors = "white",linestyles = "dashed",linewidths = 0.8)

        fig.suptitle(title)
        plt.show()
    
    return popt , pcov

# rotates axis by an angle theta
def rotation(X,Y,theta):
    V = X*np.cos(theta)+Y*np.sin(theta)
    U = -X*np.sin(theta)+Y*np.cos(theta)
    return V , U

# gaussian mathematical 1D function
def gaussian(x, A, sigma, x0):
    return 1 + np.abs(A) * np.exp(-(x-x0)**2/ ( 2*sigma**2))


# function to fit <N1> pair density
def DCE_N1(XY,x0,sigmax,Ax,offX):
    x = XY[0]-x0
    pairX = np.abs(Ax)*np.exp(-np.power(x/sigmax,2)/2) + np.abs(offX)
    return pairX

# function to fit <N2> pair density
def DCE_N2(XY,y0,sigmay,Ay,offY):
    y = XY[1]-y0
    pairY = np.abs(Ay)*np.exp(-np.power(y/sigmay,2)/2) + np.abs(offY)
    return pairY

# combination of 2D gaussian functions
def Gaussian2D(XY,x0,y0,A0,sigmax,sigmay,offset):
    x = XY[0]-x0
    y = XY[1]-y0
    f1 = np.abs(A0)*np.exp(-np.power(x/sigmax,2)/2)*np.exp(-np.power(y/sigmay,2)/2)
    return np.abs(offset) + f1

# another definition of combination of Gaussian functions
def ProdGaussian(XY,x0,y0,A0,sigmax,sigmay,A1,A2,offset):
    x = XY[0]-x0
    y = XY[1]-y0
    c0 = A0*np.exp(-np.power(x/sigmax,2)/2)*np.exp(-np.power(y/sigmay,2)/2)
    c1 = A1*np.exp(-np.power(x/sigmax,2)/2)
    c2 = A2*np.exp(-np.power(y/sigmay,2)/2)
    return c0 + c1 + c2 + np.abs(offset)


# Skewed 2D gaussian function
# see https://gregorygundersen.com/blog/2020/12/29/multivariate-skew-normal/
# see https://stackoverflow.com/questions/52975883/creating-a-multivariate-skew-normal-distribution-python
# see https://stats.stackexchange.com/questions/250874/bivariate-skewed-normal-distribution
def Skew2DGaussian(XY,x0,y0,A0,sigmaX,sigmaY,CovXY,alphax,alphay,offset):
    """ Skewed 2D gaussian function.
    Parameters
    --------------
    XY : numpy array
        axis values. XY should of the format [x,y] where x and y ate the x and y axis values.
    x0 , y0 : floats
        coordinates of the center of the distribuition
    A0 : float
        center of the distribuition
    sigmaX , sigmaY , CobXY : floats
        Covariance matrix elements
    alphax , alphay : floats
        x and y skewness values
    offset : float
        offset of the distribuition
    Returns
    --------------
        Skewed 2D gaussian function
    """
    x = XY[0]-x0
    y = XY[1]-y0
    # calculate gaussian function with covariance
    gauss = np.exp(-np.power(x/sigmaX,2)/2)*np.exp(-np.power(y/sigmaY,2)/2)*np.exp(-CovXY*x*y)
    # calculate skewness 
    skewness = 1 + erf((alphax*x + alphay*y)/np.sqrt(2))
    return A0*gauss*skewness + np.abs(offset)