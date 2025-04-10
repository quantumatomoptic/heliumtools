from heliumtools.correlations import Correlation
from heliumtools.tools import apply_ROI
from heliumtools.misc.some_plots_volume1 import heatmap_with_boxes

from scipy.optimize import curve_fit
from scipy.interpolate import RegularGridInterpolator

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
        for j, (attr_result, df_atoms, df_initial_res,axis_key, axis_integrated) in enumerate(zip(
            ["cross_1Dresult", "loc1_1Dresult", "loc2_1Dresult"],
            [self.cross, self.loc1,self.loc2], 
            [self.cross_result, self.loc1_result, self.loc2_result],
            ["Vz1+Vz2","Vz1-Vz2","Vz1-Vz2"],
            ["Vz1-Vz2","Vz1+Vz2","Vz1+Vz2"],
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

            # on dataframe with atom number per box
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
            
    def fitCorrelations(self):
        """ Function that fits local and cross correlations 
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
        fig, axes = plt.subplots(figsize = (12,9),ncols = 3, nrows=3)
           
        # for each region
        for j, (attr_result, df_atoms, df_initial_res,axis_key, axis_integrated) in enumerate(zip(
            [self.cross_1Dresult, self.loc1_1Dresult, self.loc2_1Dresult],
            [self.cross, self.loc1,self.loc2], 
            [self.cross_result, self.loc1_result, self.loc2_result],
            ["Vz1+Vz2","Vz1-Vz2","Vz1-Vz2"],
            ["Vz1-Vz2","Vz1+Vz2","Vz1+Vz2"],
            )):
            
            # check up what you do, by plotting G² and g²
            heatmap_with_boxes(ax=axes[0,j],df=df_initial_res,columns=self.var1.name,index=self.var2.name,boxes={},values="G2",cmap="Greys")
            heatmap_with_boxes(ax=axes[1,j],df=df_initial_res,columns=self.var1.name,index=self.var2.name,boxes={},values="g^2",cmap="Greys",vmax=2.2,vmin=0.9)
                
            """ Apply the selection along the integration axis """
            # on dataframe with atom number per box
            df_initial_res_roi = apply_ROI(df_initial_res, roi_list[j])

            # check up what you do, by plotting regions of integration of the 2D plots
            heatmap_with_boxes(ax=axes[0,j],df=df_initial_res_roi,columns=self.var1.name,index=self.var2.name,boxes={},values="G2",cmap=cmaps[j],cbar_bool=False,)
            heatmap_with_boxes(ax=axes[1,j],df=df_initial_res_roi,columns=self.var1.name,index=self.var2.name,boxes={},values="g^2",cmap = cmaps[j],vmax=2.2,vmin=0.9,cbar_bool=False,)
            a = df_initial_res_roi[df_initial_res_roi[axis_integrated] == df_initial_res_roi[axis_integrated].min()]
            b = df_initial_res_roi[df_initial_res_roi[axis_integrated] == df_initial_res_roi[axis_integrated].max()]
            for k in range(2):
                axes[k,j].plot(a["Vz1"],a["Vz2"],color = plt.get_cmap(cmaps[j])(1.0))
                axes[k,j].plot(b["Vz1"],b["Vz2"],color = plt.get_cmap(cmaps[j])(1.0))

           
            # check up what you do, by plotting 1D correlation
            # get rows inside the region of interest
            my_df = apply_ROI(attr_result, {axis_key: [-self.length, self.length]}).dropna()

            ax = axes[2,j]
            for jjj, Ystr in enumerate(["g^2", "g^(2) integrated","g^(2) int2"]):
                color = plt.get_cmap(cmaps[j])(0.2+0.6*jjj)
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
                    ax.plot(xth, fit_func(xth, *popt), color=color, label = "{:.2f}({:.0f})".format(maxi,100*perr))
                except Exception as e:
                    print(f"Failed to fit {e}")
                ax.grid(True, alpha = 0.5)
                ax.legend()
                ax.set_ylim([.8,2.8])
                ax.set_xlabel(axis_key)

        plt.show()
                
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

# gaussian mathematical 1D function
def gaussian(x, A, sigma, x0):
    return 1 + A * np.exp(-(x-x0)**2/ ( 2*sigma**2))

# combination of 2D gaussian functions
def Gaussian2D(XY,x0,y0,A0,sigmax,sigmay,offset,A1,A2):
    x = XY[0]
    y = XY[1]
    f1 = A0*np.exp(-np.power((x-x0)/sigmax,2))*np.exp(-np.power((y-y0)/sigmay,2))
    f2 = A1*np.exp(-np.power((x-x0)/sigmax,2))
    f3 = A2*np.exp(-np.power((y-y0)/sigmay,2))
    return offset + f1 + f2 + f3

# function to fit 2D data
def fit2D(func,df,key,xname,yname,guess,show = False):
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
    show : bool
        True if you want to plot fit results. False otherwise.
    Returns
    --------------
    popt : 1D numpy array
        optimized fit parameters found by the fit
    pcov : 2D numpy array
        covariance matrix of fit parameters
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
    zflat = zdata.flatten()
    # try to fit
    try:
        popt , pcov = curve_fit(func,xy_stack,zflat,p0=guess)
    except:
        popt = guess
        pcov = []
    
    # plot data, fit, interpolation and data-fit
    if show:
        fig, ax = plt.subplots(ncols = 4,sharey = True,sharex=True,figsize = (20,4))

        """ plot original data """
        X , Y = np.meshgrid(xdata, ydata)
        plot1 = ax[0].pcolormesh(X, Y, zdata)
        ax[0].set_xlabel('Vz1')
        ax[0].set_ylabel('Vz2')
        ax[0].set_title("data")
        cb1 = fig.colorbar(plot1,ax=ax[0])

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

        """ plot interpolation of data """

        # compute 2D interpolation
        func = RegularGridInterpolator((xdata,ydata),np.transpose(zdata),method = "cubic")
        Z = func((X1,X2))

        plot3 = ax[2].pcolormesh(X1, X2, Z)
        ax[2].set_xlabel('Vz1')
        ax[2].set_title("interpolation")
        cb3 = fig.colorbar(plot3,ax=ax[2])
        cb3.mappable.set_clim(*cb1.mappable.get_clim())

        """ plot difference between data and fit """
        Deltaz = np.abs(Gaussian2D([X,Y],*popt) - zdata)/Gaussian2D([X,Y],*popt)
        plot4 = ax[3].pcolormesh(X, Y, Deltaz)
        ax[3].set_xlabel('Vz1')
        ax[3].set_title("(data - fit)/fit")
        ax[3].set_xlim(-14,-10)
        ax[3].set_ylim(10,14)
        cb4 = fig.colorbar(plot4,ax=ax[3])

        plt.show()


    


    
    return popt , pcov