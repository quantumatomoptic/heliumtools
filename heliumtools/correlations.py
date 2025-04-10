#!/usr/bin/env python
# -*- mode:Python; coding: utf-8 -*-

# ----------------------------------
# Created on the 30-5-2022 by Victor
#
# Copyright (c) 2022 - Helium1@LCF
# ----------------------------------
#
"""
Content of correlations.py
-----------------------------

Definition of Correlation and Variable classes.
[July23] This class inherits from the DataBuilder class so that all correlation classes have the same data_builder parents.

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm, trange
import copy, random
from scipy.special import factorial
import matplotlib.cm as cm
import matplotlib.colors as colors
import time
import logging
from .data_builder import DataBuilder
from .misc.gather_data import apply_ROI, apply_ROD
from heliumtools.misc.logger import getLogger

log = getLogger(__name__)
# set the desired level of warning : DEBUG / INFO / WARNING
log.setLevel(logging.INFO)


class Correlation(DataBuilder):
    """
    Correlation class that inherits from the DatBuilder class

    Mandatory parameters
    --------------------
    atoms : pandas dataframe of 4 columns containing the cycle number, and the arrival times and positions of the atoms. The columns should be named "Cycle" ; "X" ; "Y" and "T".

    n_cycles : the number of files/cycles.

    Main attributs
    --------------------
    atoms : pandas dataframe of 4 columns containing the cycle number, and the velocities of the atoms according to the 3 axes Vx, Vy and Vz


    Other attributs
    --------------------
    boxes : dictionary giving position and size of boxes 1 and 2 for correlations.
        Forme : { "1": {"Vx": {"size": 10, "position": 0},
                        "Vy": {"size": 3, "position": 0},
                        "Vz": {"size": 0.9, "position": 56},   },
                "2": {  "Vx": {"size": 10, "position": 0},
                        "Vy": {"size": 3, "position": 0},
                        "Vz": {"size": 0.9, "position": 80},   },}

    ROI : dictionnaire, Region Of Interest
        defining a suitable ROI allows you to reduce the size of the dataframes used (in mm/s)
        WARNING: applying the apply_ROI() method acts on the atoms dataframe so if you realize that your ROI is too small, you must reload the atoms dataframe.
        Format: {"Vx": {"max":120, "min":-120}}

    ROD : dictionnaire, Region Of Desinterest
        to select only atoms outside the ROD. Same remark as for ROI and same format.

    bec_arrival_time : float, condensate arrival time in ms

    raman_kick : float, kick raman in mm/s

    var1 et var2 : Variable object (see class below), the parameters of the boxes that we are going to change to make the correlations.

    round_decimal : it turned out (LabJournal of 05/24/2022) that python does some weird rounding when it calculates Vz1 + Vz2: 
                    I therefore round all the numbers concerning Vz1 and Vz2 (or rather the self.var1.name variables) to the decimal place 
                    round_decimal (default 5) (--> see the copute_result method)


    Some Methods
    --------------------
    apply_ROI() / apply_ROD() : select atoms only in ROI / outside ROD.

    define_variable1 / define_variable2 : constructs the variable on which we will compute corrrelations. See Variable class for more infos.

    compute_correlations : construct the result dataframe in which correlations are stored.

    """

    correlation_order_max = 6

    def __init__(self, atoms, **kwargs):
        """
        Object initialization, sets parameters as the user defined, build the atoms dataframe and apply ROD and ROI.
        """
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

    def set_boxes(self, boxes):
        self.boxes = boxes.copy()

    def define_variable1(self, **kwargs):
        self.var1 = Variable(**kwargs)

    def define_variable2(self, **kwargs):
        self.var2 = Variable(**kwargs)

    def set_box_size(self, vz=0, vx=0, vy=0):
        """Useful function to set all boxes size in one line.

        Parameters
        ----------
        vz : int, optional
            size of the box along z, in mm/s, by default 0
        vx : int, optional
            size of the box along x, in mm/s, by default 0
        vy : int, optional
            size of the box along x, in mm/s, by default 0
        """
        sizes = {"Vx": vx, "Vy": vy, "Vz": vz}
        for key, value in sizes.items():
            if value:  # if the value is not zero, we set it foeach box
                for num in ["1", "2"]:
                    if "size" in self.boxes[num][key]:
                        self.boxes[num][key]["size"] = value
                    else:
                        log.warning(
                            f"[heliumtools.correlations] Setting box size failed because the box {num} format along {key} does not contain any size."
                        )

    def get_atoms_in_box(self, df, box):
        """
        Returns a dataframe with the positions of atoms inside the box.

        Parameters
        ----------
        df : atom dataframe
        box : dictionary, of the type {"Vx": {"size": 10, "position": 0}}. The dictionary entries must match the name of the 
        dataframe columns, namely Vx, Vy, Vz and Cycle. Updated since May 2023, the box dictionary entries can be 'range' 
        (array) or minimum and maximum.

        Returns
        ----------
        df : dataframe with the same columns whose atoms are all in the box.
        """
        df = apply_ROI(df, box)
        return df
        # for key, value in box.items():
        #     # Rappel : key est par ex "Vx" ; value est {"size":10, "position":0}
        #     if "range" in value:
        #         minimum = np.min(value["range"])
        #         maximum = np.max(value["range"])
        #     elif "position" in value and "size" in value:
        #         minimum = value["position"] - np.abs(value["size"]) / 2
        #         maximum = value["position"] + np.abs(value["size"]) / 2
        #     elif "minimum" in value and "maximum" in value:
        #         minimum = value["minimum"]
        #         maximum = value["maximum"]
        #     df = df[((df[key] >= minimum) & (df[key] < maximum))]

        # return df

    def merge_dataframe_on_cycles(self, df1, df2):
        """
        Merge 2 dataframe sur l'entête "Cycle". Met 0 si le df2 n'a pas de valeur à ce cycle.
        """
        df_merged = df1.merge(
            df2, how="outer", on="Cycle"
        )  # l'option "outer" permet de conserver les cycles où
        # il n'y a pas d'atomes. Pandas ajoute donc un NaN à la place.
        df_merged = df_merged.fillna(0)
        return df_merged

    def obtain_number_of_atoms_per_cycle_in_box(self, df, box, column_name="N_1"):
        """
        Retourne un dataframe de deux colonnes : une colonne avec le numéro du cycle (nom de colonne "Cycle") et une colonne avec nombre d'atomes compté dans la boite au cycle correspondant (nom de colonne column_name)

        Parameters
        ----------
            df : pandas dataframe, les données complets avec (au moins) 4 colonnes :'Cycle', 'Vx', 'Vy', 'Vz'
            box : position et taille de la boîte sur 1 à 4 axes. Exemple {"Vx": {"size": 10, "position": 0}}
            column_name : nom de la colonne du dataframe renvoyé

        Returns
        -------
            atoms_in_box : dataframe de deux colonnes avec le numéro du cycle et le nombre d'atomes dans la boîte.
            Le nom de la deuxième colonne (nombre d'atome dans la boîte) est l'argument column_name (N1 par défaut).
        """

        df = self.get_atoms_in_box(df, box)
        # We now have only atoms inside the box.
        # We count for each cycle, the number of atoms per cycle using the method
        # "value_counts". This method returns a serie and not a dataframe and I prefer retranform it
        # into a dataframe.
        atoms_in_box = (
            df.value_counts(subset="Cycle")
            .rename(column_name)
            .to_frame()
            .reset_index(inplace=False)
        )
        # atoms_in_box is now a dataframe with two columns "Cycle" and "N_1" (or column_name). However, if there were no atom at cycle 34 in the box, this cycle does not appear inside atoms_in_box. In order to have the number of atoms in the box at each cycle, we must add 0 to those cycles which does not appear.
        # cycle_dataframe is just a dataframe with n_cycles : we use it to merge and add zeros to atoms_in_box
        cycle_dataframe = pd.DataFrame(self.cycles_array, columns=["Cycle"])
        atoms_in_box = self.merge_dataframe_on_cycles(cycle_dataframe, atoms_in_box)
        return atoms_in_box

    def counts_atoms_in_boxes_one_variable(self, df, var, box, column_name = "N_1"):
        """
        Takes as arguments a dataframe of atoms, a variable and a box. 
        For each "value" of "variable", it redefines the size/position of the box and retrieves the number of atoms in the box each cycle. 
        It returns a dataframe of 3 columns: one with the cycles named "Column", 
        one with the number of atoms in the given cycle (named column_name) 
        and one with the value of the position/size of the box (named var .name)

        Parameters
        ----------
        df : pandas dataframe, dataframe with atoms
        var : Variable
        box : dictionary, position and size of the box on 1 to 4 axes.
        Example {"Vx": {"size": 10, "position": 0}}
        column

        Returns
        ----------
        result : pandas dataframe
            dataframe with 3 columns:
                "Cycle": with the cycle number
                column_name (ex: "N1"): with the number of atoms in the box
                var.name (ex: "ΔVx"): with the value of the position/size of the box.
        """
        # We go through the different values ​​of the Variable var (tqdm --> waiting bar)
        # for i in tqdm(range(var.n_step), desc="Gathering {}".format(var.name)):
        for i in range(var.n_step):
            # We change the box according to the i-th value of var
            box[var.axe][var.type] = var.get_value_i(i)
            # We retrieve the dataframe with the number of atoms in the box at each cycle
            dataframe = self.obtain_number_of_atoms_per_cycle_in_box(
                df, box, column_name=column_name
            )
            dataframe[var.name] = var.get_value_i(i)
            # We add the column with the variable value
            if i == 0:
                result = dataframe
            else:
                result = pd.concat([result, dataframe])
        # the indexes repeat themselves: it's meh: I reset them.
        result.reset_index(drop=True)
        return result

    ###########################################################
    ######## MANAGEMENT OF CASES BEFORE CALCULATION  #######
    ###########################################################
    def compute_correlations(self):
        """
        
        This function computes the correlations and manages the different scanning scenarios.        
        
        """

        # Case 1: we do not scan any parameters, we just want the correlation between two boxes.
        if (self.var1 == None) and (self.var2 == None):
            atoms_box1 = self.obtain_number_of_atoms_per_cycle_in_box(
                self.atoms, self.boxes["1"], column_name="N_1"
            )
            atoms_box2 = self.obtain_number_of_atoms_per_cycle_in_box(
                self.atoms, self.boxes["2"], column_name="N_2"
            )
            total_atoms = self.merge_dataframe_on_cycles(atoms_box1, atoms_box2)
            self.total = total_atoms
            corr_names = self.quantity_of_interest()
            corr_values = self.quantity_of_interest(total_atoms)
            self.result = pd.DataFrame([corr_values], columns=corr_names)

        # Case 2: we scan only one parameter: in this case, we define the second scanned variable with a single parameter 
        # by calling define_new_variable_with_one_value then we call the compute_correlations method again to be in case 3
        elif (self.var1 == None) and (self.var2 != None):
            # in this case we will define the variable 1
            # print("I defined myself variable1")
            self.define_new_variable_with_one_value(self.var2, 1)
            self.compute_correlations()

        elif (self.var1 != None) and (self.var2 == None):
            # print("I defined myself variable2")
            self.define_new_variable_with_one_value(self.var1, 2)
            self.compute_correlations()

        # Case 3: we scan parameters belonging to two different boxes
        elif self.var1.box != self.var2.box:
            self.compute_correlations_different_box_scanned()

        # Case 4: we scan parameters belonging to the same box
        else:
            self.compute_correlations_same_box_scanned()

    def compute_correlations_same_box_scanned(self):
        """
        This method is called when the two variables scan the same ax. Built a total dataframe and call the method compute_result.
        """
        boxes = ["1", "2"]
        if self.var1.box not in boxes:
            print(
                f"[ERROR] : the var1 box is neither '1' nor '2' but {self.var1.box } "
            )
            return
        boxes.remove(self.var1.box)
        box_not_scanned = boxes[0]
        # The following is a dataframe with only two column wich are the number of atom in the box that is not scanned
        atom_beam_not_scanned = self.obtain_number_of_atoms_per_cycle_in_box(
            self.atoms, self.boxes[box_not_scanned], column_name="N_" + box_not_scanned
        )
        result_var_list = []
        for i in range(self.var1.n_step):
            # We must set the new value for the box which is scanned :
            self.boxes[self.var1.box][self.var1.axe][self.var1.type] = (
                self.var1.get_value_i(i)
            )
            df = self.counts_atoms_in_boxes_one_variable(
                self.atoms,
                self.var2,
                self.boxes[self.var2.box],
                column_name="N_" + self.var1.box,
            )
            df[self.var1.name] = self.var1.get_value_i(i)
            result_var_list.append(df)
        result_var = pd.concat(result_var_list).reset_index()
        # print(len(result_var_2))
        # result_var2_list = [result_var_2]

        self.total = pd.merge(atom_beam_not_scanned, result_var, on="Cycle")
        self.compute_result(self.total)

    def compute_correlations_different_box_scanned(self):
        """
        Method for calculating correlations when var1 and var2 (the scanned parameters) correspond to two different boxes.
        """

        """ To understand the algorithm consider the example: we want the correlation map in the Vz direction for both pairs, that is 
        determine g^(2)(Vz1 , Vz2). In this case, the scanned variable is Vz, that is we create several boxes in Vz, for both beams in the
        pairs, and cound the number of atoms in each box. Notice however, that the scanned variable can be Vx, for example, and they do not
        need to be the same for both beams, e.g, for va1 the scanned variable can be Vz and for var2 the scanned varibale can be Vx. """

        """ STEP 1: we retrieve the number of atoms in the boxes associated with var1 and var2. """

        # Let's say that intuitively, var1 corresponds to box 1 and var2 to box 2 but this is not necessary in the code.

        # We start with var1. We start by retriving all atoms present in the box transverse to the scanned direction.
        # E.g if Vz is the scanned direction, then we collect all atoms in the transverse (Vx,Vy) box defined in var1 
        # This lightens the calculation.

        # create copy of box defined in var1
        box = self.boxes[self.var1.box].copy()
        # remove the scanned axis of the box, e.g Vz
        posi_and_size = box.pop(self.var1.axe)
        # create scanned box
        scanned_box = {self.var1.axe: posi_and_size}

        # Now we have two boxes: 
        # scanned_box -> corresponds to the scanned direction var1.axe, e.g Vz. 
        # box -> transverse box, along transverse directions, e.g Vx and Vy. This box is the same for all different positions of var1.axe.

        # get all atoms that are inside the transverse box, e.g Vx and Vy.
        df_atoms_var1 = self.get_atoms_in_box(self.atoms, box)
        
        # we scan the scanned direction: we split it into several small boxes and count the number of atoms in each box.
        # the size of each box is determined by var1
        # The return is result_var1 which is a dataframe with the number of atoms in each small box at each cycle. 
        # result_var1 has 3 columns with headers "Cycle", "N_i" with i = 1 or 2 depending on the box and the name of the scanned variable e.g. "ΔVz". 
        # See the documentation of counts_atoms_in_boxes_one_variable for more details.
        result_var1 = self.counts_atoms_in_boxes_one_variable(df_atoms_var1,self.var1, scanned_box,column_name="N_"+self.var1.box)
        
        ## check 21 of may:
        column_name = "N_" + self.var1.box
        df = result_var1[result_var1[column_name] < 0]
        if len(df) > 1:
            print(df)

        # We do the same with var2
        # create copy of box defined in var2
        box = self.boxes[self.var2.box].copy()
        # remove the scanned axis of the box
        posi_and_size = box.pop(self.var2.axe)
        # create scanned box
        scanned_box = {self.var2.axe: posi_and_size}
        # get all atoms that are inside the transverse box
        df_atoms_var2 = self.get_atoms_in_box(self.atoms, box)
        # split scanned direction into several small boxes and count the number of atoms in each box.
        result_var2 = self.counts_atoms_in_boxes_one_variable(df_atoms_var2, self.var2, scanned_box, column_name="N_" + self.var2.box)

        column_name = "N_" + self.var2.box
        df = result_var2[result_var2[column_name] < 0]
        if len(df) > 1:
            print(df)
        
        """ STEP2: merge dataframes """

        # We build the total dataframe, which initially contains 5 columns: Cycle, N_1 and N_2 (the number of atoms in box 1 and 2) and self.var1 and self.var2 (the position/size of the boxes during the scan).
        # The number of rows in total is therefore (Number_of_cycles) * (Number_of_different_var1) * (Number_of_different_var2).

        total = pd.merge(result_var1, result_var2, on = "Cycle")

        """ STEP3: computes quantities of interest """
        self.compute_result(total)

    def compute_opposite_momenta_correlations(self):
        """
        Compute correlations in self.result when you only want to scan one variable and make correlations between k and -k. 
        The parameters that are moved must be defined in var1.
        """

        # set var1 as master
        box = self.var1.box
        
        if box == "1":
            box = "2"
        else:
            box = "1"
        
        # define var2 so that its box lies in the opposite momentum region
        axe = self.var1.axe
        type = self.var1.type
        name = self.var1.name
        values = self.var1.values
        self.define_variable2(box=box, axe=axe, type=type, name="-" + name, values=-1 * values)

        # Now we do the same thing as in compute_correlations_different_box_scanned()

        """ STEP 1: we retrieve the number of atoms in the boxes associated with var1 and var2. 
        Let's say that intuitively, var1 corresponds to box 1 and var2 to box 2 but this is not necessary in the code. """

        # Start with var1
        # We only retrieve the atoms present in the box according to the two axes that do not vary. 
        # For example, if we are varying the position of the box according to Vx, we retrieve the atoms that already verify the right conditions according to Vy and Vz to lighten the calculations.
        box = self.boxes[self.var1.box].copy()
        # We remove the axis concerned by the scan
        posi_and_size = box.pop(self.var1.axe)
        # the axis concerned by the scan is therefore
        scanned_box = {self.var1.axe: posi_and_size}
        # So we have two boxes now: one with two axes (box) that are not modified at each box position/size and another (scanned_box) with a single axis that corresponds to the scanned axis var1.axis
        df_atoms_var1 = self.get_atoms_in_box(self.atoms, box)
        # Result var1 is a dataframe with the number of atoms in the box at each cycle for the different positions of the box. 
        # So it has 3 columns with headers "Cycle", "N_i" with i = 1 or 2 depending on the box and the name of the scanned variable for example "ΔVx". 
        # See the documentation of counts_atoms_in_boxes_one_variable for more details.
        result_var1 = self.counts_atoms_in_boxes_one_variable(df_atoms_var1, self.var1, scanned_box, column_name="N_" + self.var1.box)

        # We do the same with var2
        box = self.boxes[self.var2.box].copy()
        posi_and_size = box.pop(self.var2.axe)
        scanned_box = {self.var2.axe: posi_and_size}
        df_atoms_var2 = self.get_atoms_in_box(self.atoms, box)
        result_var2 = self.counts_atoms_in_boxes_one_variable(df_atoms_var2, self.var2, scanned_box, column_name="N_" + self.var2.box)

        """  STEP2  """
        result_var1["abs(varname)"] = np.abs(result_var1[self.var1.name])
        result_var2["abs(varname)"] = np.abs(result_var2[self.var2.name])

        # We now merge the dataframes on two keys and no longer on just one.
        total = pd.merge(result_var1, result_var2)

        """ STEP3 : computes quantites of interest """
        self.compute_result(total) 

    def define_new_variable_with_one_value(self, var, new_var_number):
        """ This function generates variable 1 or 2 (corresponding to var_number) by copying the parameters of var. 
        For example, if we scan variable 1, the size of boxes according to Vz, this function will define variable 2 as scanning the size of 
        the boxes in the box opposite to variable 1 with a single parameter, the one defined in box. 
        This allows to systematically use the compute_correlations_different_box_scanned function even when we have only one box defined.

        Parameters
        ----------
        var : object of the Variable class

        var_number : integer 1 or 2
            If it is 1, we will define a new variable 1, if it is 2, we will define the new variable 2.
        """
        if new_var_number not in [1, 2]:
            raise (
                "Error : this variable number is not defined ({}). Choose either 1 or 2.".format(
                    new_var_number
                )
            )
        # Pour définir une variable, il faut
        # --> le numéro de sa boîte : "1"  ou "2". C'est celui opposé à var.
        liste = ["1", "2"]
        liste.remove(var.box)
        box_number = liste[0]
        # --> son axe et son type, qui sont le même que la boîte déjà définie.
        axe = var.axe  # l'axe du paramètre scanné (Vx, Vy ou Vz)
        type = var.type  # le type : size ou position
        # On récupère son nom en prenant les paramètres par défaut.
        name = var.built_name(add_box_number=False) + box_number
        # si jamais le nom de la variable correspond à l'autre, cela va poser problème donc on rajoute un _default.
        if name == var.name:
            name += "_default"
        values = [self.boxes[box_number][axe][type]]
        if new_var_number == 1:
            self.define_variable1(
                box=box_number, axe=axe, type=type, name=name, values=values
            )
        elif new_var_number == 2:
            self.define_variable2(
                box=box_number, axe=axe, type=type, name=name, values=values
            )

    def quantity_of_interest(self, dataframe=""):
        """
        Takes as argument a dataframe whose columns are the number of atoms in boxes 1 and 2 at each cycle1 
        Returns a list of floats with the different values ​​of interest. If dataframe = "", returns a list of strings with the name of each calculated quantity.
        To add values ​of interest to calculate, simply add the name of the variable to the column_names list and add it to result.

        Parameters
        ----------
        dataframe : pandas dataframe with (at least) 3 columns
        "Cycle" the cycle number
        "N_1", with the number of atoms in box 1
        "N_2", with the number of atoms in box 2

        Returns
        ----------
        column_names (if no dataframe is given): list of strings with the name of the calculated values
        computed_quantities (if a dataframe is given): list of floats with the calculated quantity of interest.
        """
        
        column_names = [
            "N_1",
            "N_2",
            "N_1+N_2",
            "variance",
            "normalized variance",
            "g^2",
        ]

        if type(dataframe) == str:
            return column_names

        N1 = np.sum(dataframe["N_1"]) / self.n_cycles
        N2 = np.sum(dataframe["N_2"]) / self.n_cycles

        variance = (
            np.sum((dataframe["N_1"] - dataframe["N_2"]) ** 2) / self.n_cycles
            - (np.sum(dataframe["N_1"] - dataframe["N_2"]) / self.n_cycles) ** 2
        )
        normalized_variance = variance / (N1 + N2)

        # Calcul de g2
        numerator = np.sum(dataframe["N_1"] * dataframe["N_2"]) / self.n_cycles
        denominator = (
            np.sum(dataframe["N_1"]) * np.sum(dataframe["N_2"]) / (self.n_cycles**2)
        )
        if denominator > 0:
            g2 = numerator / denominator
        else:
            g2 = 1

        computed_quantities = [
            N1,
            N2,
            N1 + N2,
            variance,
            normalized_variance,
            g2,
        ]
        return computed_quantities

    ###########################################################
    ############### FONCTION DE CALCUL FINALE #################
    ###########################################################

    def compute_result(self, total):
        """
        Compute the correlation dataframe result using total. 
        Total is a dataframe with for each cycle and each box position/size the number of atoms in box 1 and in box 2.

        /!\ /!\ Only works if we have two variables --> to do for a 1D scan? or we keep the naive way because it is not long.

        Parameters
        ----------
        total : pandas dataframe of 5 columns: 'Cycle' the cycle, 'N_1' and 'N_2' number of atoms in box 1 and 2, self.var1 and self.var2 the position/size of the boxes when scanned. 
        The number of rows in total is therefore Number_of_cycles x Number_of_different_var1 x Number_of_different_var2.

        Built
        ----------
        self.result : pandas dataframe of Number_of_different_var1 x Number_of_different_var2 rows. Its different columns are
        var1.name and var2.name : the name of the scanned variables
            "N_1", "N_2" :  average over cycles of the number of atoms in box 1, 2
            "N_1+N_2" : sum of average populations
            "N_1-N_2" : mean population difference
            "(N_1-N_2)^2" : average (over cycles) of the squared deviations
            "N_1*N_2" : average of the product of the number of atoms
            "variance" : variance, this is the quantity <(N1 - N2)^2> / <N1-N2>^2 where the average is taken over the cycles.
            "normalized variance" : normalized variance, variance divided by the sum of the populations: variance / <N1 + N2>
            "<N_1>*<N_2>" : product of the average number of atoms
        """
        
        """ compute some usefull quantities for future calculation and add to dataset"""
        total["N_1*N_2"] = total["N_1"] * total["N_2"]
        total["N_1**2"] = total["N_1"] ** 2
        total[":N_1**2:"] = total["N_1"] ** 2 - total["N_1"]
        total["N_2**2"] = total["N_2"] ** 2
        total[":N_2**2:"] = total["N_2"] ** 2 - total["N_2"]
        total[":N_1^2xN_2^2:"] = total["N_1"] * (total["N_1"] - 1) * total["N_2"] * (total["N_2"] - 1)
        total["N_1-N_2"] = total["N_1"] - total["N_2"]
        total["(N_1-N_2)^2"] = (total["N_1"] - total["N_2"]) ** 2
        total["N_1+N_2"] = total["N_1"] + total["N_2"]
        total["(N_1-N_2)^4"] = (total["N_1"] - total["N_2"]) ** 4
        total["M^2 jasukula"] = total["(N_1-N_2)^2"] / (total["N_1+N_2"])
        total["M jaskula"] = total["N_1-N_2"] / np.sqrt(total["N_1+N_2"])
        total["[N_1*N_2]**2"] = total["N_1*N_2"] * total["N_1*N_2"]
        total[":N_1**2:**2"] = total[":N_1**2:"] ** 2
        total[":N_2**2:**2"] = total[":N_2**2:"] ** 2
        total["N_1**2*N_2**2"] = total["N_1"] ** 2 * total["N_2"] ** 2
        total["N_1**2*N_2"] = total["N_1"] ** 2 * total["N_2"]
        total["N_1*N_2**2"] = total["N_1"] * total["N_2"] ** 2
        total["Jz fluctu"] = 1/2*(np.sqrt((total["N_1+N_2"] - 1) * np.heaviside(total["N_1+N_2"], 0))* total["N_1-N_2"])
        total["Jz fluctu^2"] = total["Jz fluctu"] ** 2

        """ compute nth order correlations """
        # 1th order correlation
        total[":N_1^1:"] = total["N_1"]
        total[":N_2^1:"] = total["N_2"]
        # nth order correlation
        for i in range(2, self.correlation_order_max + 1):
            total[f":N_1^{i}:"] = total[f":N_1^{i-1}:"] * (total["N_1"] - i + 1)
            total[f":N_2^{i}:"] = total[f":N_2^{i-1}:"] * (total["N_2"] - i + 1)

        # update object dataset
        self.total = total
        
        # we average the data over the cycles (we therefore group them by different values ​​of var1 and var2)
        self.result = total.groupby([self.var1.name, self.var2.name], as_index=False).mean()
        # compute std
        if self.compute_errors:
            self.errors = total.groupby([self.var1.name, self.var2.name], as_index=False).std()
        column_name = "N_1*N_2"
        df = self.result[self.result[column_name] < 0]
        if len(df) > 1:
            print("there is an issue !!!!!")
            print(df)

        # # On fait ensuite une petite manipulation pour calculer l'erreur sur la variance.
        # # L'idée est de définir la variable (N_1-N_2)^2-moy(N_1-N_2)^2 puis de dire à l'ordinateur de calculer sa variance tout seul pour éviter de mettre la formule très longue et compliquée (wiki du 2 juin 2022). On reconstruit un dataframe avec les numéros de cycles
        # df1 = self.result[[self.var1.name, self.var2.name, "N_1-N_2"]]
        # df2 = pd.DataFrame({"Cycle": np.linspace(1, self.n_cycles, self.n_cycles)})
        # new_df = pd.merge(df1, df2, how="cross")
        # new_df = new_df[["Cycle", self.var1.name, self.var2.name, "N_1-N_2"]]
        # new_df.columns = new_df.columns.str.replace("N_1-N_2", "mean(N_1-N_2)")
        # # new_df est donc un dataframe avec 4 colonnes : une kz1, une kz2, une avec le cycle et une avec la moyenne (N1-N2). Bien entendu, à chaque cycle la moyenne N1-N2 est la même. NB : kz1 est de façon générale var1.name mais souvent kz1.
        # # On veut ajouter au dataframe total la colonne moy(N1-N2).
        # total = pd.merge(total, new_df)
        # total["(N_1-N_2)^2-mean(N_1-N_2)^2"] = (
        #     total["(N_1-N_2)^2"] - total["mean(N_1-N_2)"] ** 2
        # )
        # if self.compute_errors:
        #     error = total.groupby(
        #         [self.var1.name, self.var2.name], as_index=False
        #     ).std()

        # print("Total dataframe is summed already")

        """ compute Variance, g^2 , ... """

        # ---------------
        # Variance
        # ---------------
        self.result["variance"] = self.result["(N_1-N_2)^2"] - self.result["N_1-N_2"]**2
        self.result["normalized variance"] = self.result["variance"] / (self.result["N_1"] + self.result["N_2"])

        # ---------------
        # Calculate de Denis's criterion (no projector aka no post-selection)
        # ---------------
        self.result["denis2"] = (self.result["Jz fluctu^2"] - self.result["Jz fluctu"] ** 2) / self.result["N_1*N_2"]

        # ---------------
        # Calculate de g^2
        # ---------------
        self.result["g^2"] = self.result["N_1*N_2"] / (self.result["N_1"] * self.result["N_2"])
        
        # -----------------
        # Calculate de g^n local
        # -----------------
        for i in range(2, self.correlation_order_max + 1):
            for j in [1, 2]:
                self.result[f"g_{j}^({i})"] = self.result[f":N_{j}^{i}:"] / self.result[f"N_{j}"] ** i
                if self.compute_errors:
                    self.result[f"U(g_{j}^({i}))"] = (
                        self.result[f"g_{j}^({i})"]
                        * np.sqrt(
                            self.errors[f":N_{j}^{i}:"] ** 2
                            / self.result[f":N_{j}^{i}:"] ** 2
                            + i
                            * self.errors[f"N_{j}"] ** 2
                            / self.result[f"N_{j}"] ** 2
                        )
                        / np.sqrt(self.n_cycles)
                    )
                else:
                    ## here we miss out the standard deviation of N1**i
                    self.result[f"U(g_{j}^({i}))"] = (
                        self.result[f"g_{j}^({i})"]
                        * np.sqrt(
                           i* (self.result[f"N_{j}"] ** 2 + self.result[f"N_{j}"])
                            / self.result[f"N_{j}"] ** 2
                        )
                        / np.sqrt(self.n_cycles)
                    )

        # -------------------------------
        # Calculs de delta, qty Parentani
        # -------------------------------
        self.result["Delta"] = - self.result["N_1*N_2"] / 2 + self.result["N_1"] * self.result["N_2"]
        self.result["-Delta"] = - self.result["Delta"]

        # ---------------
        # Calculs de g^4
        # ---------------
        self.result["g^4"] = self.result[":N_1^2xN_2^2:"] / (self.result["N_1"] ** 2 * self.result["N_2"] ** 2)

        # normalement ça doit être les mêmes...
        self.result["g^4 bis"] = (
            self.result["N_1**2*N_2**2"]
            - self.result["N_1*N_2**2"]
            - self.result["N_1**2*N_2"]
            + self.result["N_1*N_2"]
        ) / (self.result["N_1"] ** 2 * self.result["N_2"] ** 2)
        ## define the minimum of the fourth order correlation function for thermal gaussian state
        self.result["g^4 mini"] = (
            16 * self.result["g^2"] + 4 * (self.result["g^2"] - 1) ** 2 - 12
        )
        self.result["g^4 maxi"] = (
            16 * self.result["g^2"] + 6 * (self.result["g^2"] - 1) ** 2 - 12
        )
        g2 = self.result["g^2"]
        g4 = self.result["g^4"]
        # Calcul de theta
        self.result["theta_g^4"] = (g4 - (16 * g2 + 4 * (g2 - 1)**2 - 12)) / (2 * (g2 - 1)**2)
         
        # ---------------
        # Calculs de corrélations locales
        # ---------------
        self.result[":N_1**2:"] = self.result["N_1**2"] - self.result["N_1"]
        self.result[":N_2**2:"] = self.result["N_2**2"] - self.result["N_2"]

        # ---------------
        # Cuachy-Schwarz
        # ---------------
        self.result["C-S"] = self.result["N_1*N_2"] / np.sqrt((self.result["N_1**2"] - self.result["N_1"])* (self.result["N_2**2"] - self.result["N_2"]))

        self.result["C-S difference"] = self.result["N_1*N_2"] - np.sqrt((self.result["N_1**2"] - self.result["N_1"])* (self.result["N_2**2"] - self.result["N_2"]))

        self.result["G^2(k1,k1)"] = self.result["N_1**2"] - self.result["N_1"]
        self.result["G^2(k2,k2)"] = self.result["N_2**2"] - self.result["N_2"]
        self.result["g^2(k1,k1)"] = (
            self.result["N_1**2"] - self.result["N_1"]
        ) / self.result["N_1"] ** 2
        self.result["g^2(k2,k2)"] = (
            self.result["N_2**2"] - self.result["N_2"]
        ) / self.result["N_2"] ** 2
        self.result["var(M jaskula)"] = (
            self.result["M^2 jasukula"] - self.result["M jaskula"] ** 2
        )

        # we remove the shot noise if requested by the user.
        if self.remove_shot_noise:
            if self.var1.type == "position":
                VJ = self.var1.axe  # I assume that we scanned the same axe
                local_condition = self.result[self.var1.name] == self.result[self.var2.name]
                not_scanned_axes = ["Vx", "Vy", "Vz"]
                not_scanned_axes.remove(self.var1.axe)
                if self.var1.axe != self.var2.axe:
                    print(
                        "[WARNING] THIS IS NOT YET TAKE INTO ACCOUNT IN THE CODE. PLEASE CHANGE ME."
                    )
                elif (
                    self.boxes["1"][not_scanned_axes[0]]["position"]
                    != self.boxes["2"][not_scanned_axes[0]]["position"]
                ) or (
                    self.boxes["1"][not_scanned_axes[1]]["position"]
                    != self.boxes["2"][not_scanned_axes[1]]["position"]
                ):
                    # if len(local_condition) > 0:
                    #     print(
                    #         "[WARNING] Shot Noise has not been taken off weird because boxes are do not have the same center. "
                    #     )
                    pass
                elif (
                    (self.boxes["1"]["Vz"]["size"] != self.boxes["2"]["Vz"]["size"])
                    or (self.boxes["1"]["Vy"]["size"] != self.boxes["2"]["Vy"]["size"])
                    or (self.boxes["1"]["Vx"]["size"] != self.boxes["2"]["Vx"]["size"])
                ):
                    pass
                    # if len(local_condition) > 0:
                    #     print(
                    #         "[WARNING] Shot Noise has not been taken off weird because boxes do not have the same size. Please be carefull when delaing with local correlations !"
                    #     )
                else:
                    self.result.loc[local_condition, "g^2"] = (self.result["N_1*N_2"] - self.result["N_1"]) / (self.result["N_1"] * self.result["N_2"])
                    self.result.loc[local_condition, "N_1*N_2"] = self.result["N_1*N_2"] - self.result["N_1"]
                    ## Recompute Cauchy-Schwarz
                    self.result["C-S"] = self.result["N_1*N_2"] / (
                        np.sqrt(
                            (self.result["N_1**2"] - self.result["N_1"])
                            * (self.result["N_2**2"] - self.result["N_2"])
                        )
                    )

                    self.result["C-S difference"] = self.result["N_1*N_2"] - (
                        np.sqrt(
                            (self.result["N_1**2"] - self.result["N_1"])
                            * (self.result["N_2**2"] - self.result["N_2"])
                        )
                    )
            if self.var1.type == "size":
                if (
                    (
                        self.boxes["1"]["Vz"]["position"]
                        == self.boxes["2"]["Vz"]["position"]
                    )
                    and (
                        self.boxes["1"]["Vy"]["position"]
                        == self.boxes["2"]["Vy"]["position"]
                    )
                    and (
                        self.boxes["1"]["Vx"]["position"]
                        == self.boxes["2"]["Vx"]["position"]
                    )
                ):
                    # le shot noise correspond à la plus petite valeur entre le nombre d'atome dans la boite 1 et le nombre d'atomes dans la boite 2.
                    mini_N1_N2 = 0.5 * (
                        self.result["N_1"]
                        + self.result["N_2"]
                        - np.abs(self.result["N_1"] - self.result["N_2"])
                    )
                    self.result["N_1*N_2"] = self.result["N_1*N_2"] - mini_N1_N2
                    self.result["g^2"] = self.result["N_1*N_2"] / (
                        self.result["N_1"] * self.result["N_2"]
                    )

        """ Compute standard deviations of g2, variance, ... """

        self.result["N_1 std"] = np.sqrt(
            self.result["N_1**2"] - self.result["N_1"] ** 2
        )
        self.result["N_2 std"] = np.sqrt(
            self.result["N_2**2"] - self.result["N_2"] ** 2
        )
        self.result["N_1 rel"] = self.result["N_1 std"] / self.result["N_1"]
        self.result["N_2 rel"] = self.result["N_2 std"] / self.result["N_2"]
        self.result["N_1*N_2 std"] = np.sqrt(
            self.result["[N_1*N_2]**2"] - self.result["N_1*N_2"] ** 2
        ) / np.sqrt(self.n_cycles)

        self.result["N_1-N_2 std"] = self.result["N_1-N_2"] * np.sqrt(
            self.result["N_1 rel"] ** 2 + self.result["N_2 rel"] ** 2
        )
        self.result["(N_1-N_2)^2 std"] = self.result["(N_1-N_2)^2"] * np.sqrt(
            self.result["N_1 rel"] ** 2 + self.result["N_2 rel"] ** 2
        )
        self.result["N_1+N_2 std"] = self.result["(N_1-N_2)^2"] * np.sqrt(
            self.result["N_1 rel"] ** 2 + self.result["N_2 rel"] ** 2
        )
        self.result["variance std"] = self.result["variance"] * np.sqrt(
            self.result["N_1 rel"] ** 2 + self.result["N_2 rel"] ** 2
        )
        #### Defining error
        self.result["N_1 error"] = self.result["N_1 std"] / np.sqrt(self.n_cycles)
        self.result["N_2 error"] = self.result["N_2 std"] / np.sqrt(self.n_cycles)

        # g² and cauchy schwarz
        self.result["N_1*N_2 error"] = self.result["N_1*N_2 std"] / np.sqrt(
            self.n_cycles
        )

        self.result["g^2 error"] = np.sqrt(
            (self.result["N_1 error"] / self.result["N_1"]) ** 2
            + (self.result["N_2 error"] / self.result["N_2"]) ** 2
            + (self.result["N_1*N_2 error"] / self.result["N_1"] / self.result["N_2"])
            ** 2
        )
        self.result[":N_1**2: error"] = np.sqrt(
            (self.result[":N_1**2:**2"] - self.result[":N_1**2:"] ** 2) / self.n_cycles
        )
        self.result[":N_2**2: error"] = np.sqrt(
            (self.result[":N_2**2:**2"] - self.result[":N_2**2:"] ** 2) / self.n_cycles
        )
        self.result["C-S error"] = np.sqrt(
            (self.result["N_1*N_2 error"] / self.result["N_1"] / self.result["N_2"])
            ** 2
            + (self.result[":N_1**2: error"] / self.result["N_1"] ** 2) ** 2
            + (self.result[":N_2**2: error"] / self.result["N_2"] ** 2) ** 2
        )
        self.result["C-S difference error"] = np.sqrt(
            (self.result["N_1*N_2 error"]) ** 2
            + (self.result[":N_1**2: error"]) ** 2
            + (self.result[":N_2**2: error"]) ** 2
        )

        ## Variance
        self.result["N_1-N_2 error"] = self.result["N_1-N_2 std"] / np.sqrt(
            self.n_cycles
        )
        self.result["(N_1-N_2)^2 error"] = self.result["(N_1-N_2)^2 std"] / np.sqrt(
            self.n_cycles
        )
        self.result["N_1+N_2 error"] = self.result["N_1+N_2 std"] / np.sqrt(
            self.n_cycles
        )
        self.result["variance error"] = self.result["variance std"] / np.sqrt(
            self.n_cycles
        )

        self.result["normalized variance error"] = self.result[
            "normalized variance"
        ] * np.sqrt(
            (self.result["variance error"] / self.result["variance"]) ** 2
            + (self.result["N_1+N_2 error"] / self.result["N_1+N_2"]) ** 2
        )

        self.result["g^2 error"] = self.result["g^2"] * np.sqrt(
            (self.result["N_1*N_2 error"] / self.result["N_1*N_2"]) ** 2
            + (self.result["N_1 error"] / self.result["N_1"]) ** 2
            + (self.result["N_2 error"] / self.result["N_2"]) ** 2
        )
        if self.remove_shot_noise:
            local_condition = self.result[self.var1.name] == self.result[self.var2.name]
            # il y a certainement une erreur ici --> modifier pour prendre en compte le shot noise justement
            self.result["g^2 error"] = self.result["g^2"] * np.sqrt(
                (self.result["N_1*N_2 error"] / self.result["N_1*N_2"]) ** 2
                + 2 * (self.result["N_1 error"] / self.result["N_1"]) ** 2
                + (self.result["N_2 error"] / self.result["N_2"]) ** 2
            )

        self.result["normalized variance error2"] = (
            (self.result["(N_1-N_2)^4"] - self.result["(N_1-N_2)^2"] ** 2)
            / self.n_cycles
            / self.result["N_1+N_2"]
        )

        # ---------------
        # Densities are expressed with units
        # ---------------
        try:
            self.result["N_1 (at/(mm/s)^3)"] = self.result["N_1"] / (
                self.boxes["1"]["Vx"]["size"]
                * self.boxes["1"]["Vy"]["size"]
                * self.boxes["1"]["Vz"]["size"]
            )
            self.result["N_2 (at/(mm/s)^3)"] = self.result["N_2"] / (
                self.boxes["2"]["Vx"]["size"]
                * self.boxes["2"]["Vy"]["size"]
                * self.boxes["2"]["Vz"]["size"]
            )
            self.result["N_1 (at/(mm/s)^3) std"] = self.result["N_1 std"] / (
                self.boxes["1"]["Vx"]["size"]
                * self.boxes["1"]["Vy"]["size"]
                * self.boxes["1"]["Vz"]["size"]
            )
            self.result["N_2 (at/(mm/s)^3) std"] = self.result["N_2 std"] / (
                self.boxes["2"]["Vx"]["size"]
                * self.boxes["2"]["Vy"]["size"]
                * self.boxes["2"]["Vz"]["size"]
            )
            self.result["N_1 (at/(mm/s)^3) error"] = self.result["N_1 error"] / (
                self.boxes["1"]["Vx"]["size"]
                * self.boxes["1"]["Vy"]["size"]
                * self.boxes["1"]["Vz"]["size"]
            )
            self.result["N_2 (at/(mm/s)^3) error"] = self.result["N_2 error"] / (
                self.boxes["2"]["Vx"]["size"]
                * self.boxes["2"]["Vy"]["size"]
                * self.boxes["2"]["Vz"]["size"]
            )
        except KeyError:
            pass
        # ---------------
        # We add the difference and the average of var1 and var2
        # ---------------
        self.result[f"({self.var1.name}+{self.var2.name})/2"] = np.round(
            (self.result[self.var1.name] + self.result[self.var2.name]) / 2,
            self.round_decimal,
        )
        self.result[f"{self.var1.name}+{self.var2.name}"] = np.round(
            (self.result[self.var1.name] + self.result[self.var2.name]),
            self.round_decimal,
        )
        self.result[f"({self.var1.name}-{self.var2.name})/2"] = np.round(
            (self.result[self.var1.name] - self.result[self.var2.name]) / 2,
            self.round_decimal,
        )
        self.result[f"{self.var1.name}-{self.var2.name}"] = np.round(
            (self.result[self.var1.name] - self.result[self.var2.name]),
            self.round_decimal,
        )
        self.result[f"({self.var2.name}-{self.var1.name})/2"] = np.round(
            (self.result[self.var2.name] - self.result[self.var1.name]) / 2,
            self.round_decimal,
        )
        self.result[f"{self.var2.name}-{self.var1.name}"] = np.round(
            (self.result[self.var2.name] - self.result[self.var1.name]),
            self.round_decimal,
        )
        try:
            self.add_denis_criterion(self.total)
        except Exception as e:
            pass
            # log.error(
            #     f"[Correlations] Failed to add Denis's quantum criteria. Error is : {e}"
            # )

    def add_denis_criterion(self, total):
        #### Denis's criteria
        # keep only a small number of columns for memory usage
        to_keep = total[
            [self.var1.name, self.var2.name, "N_1+N_2", "N_1-N_2", "N_1*N_2"]
        ]
        # choose only shots for which the number of detected atoms is greater than 0
        post_selec = copy.deepcopy(to_keep[to_keep["N_1+N_2"] > 0])
        # define sqrt(N-1)*Jz
        post_selec["Jz fluctu"] = (
            1 / 2 * (np.sqrt(post_selec["N_1+N_2"] - 1) * post_selec["N_1-N_2"])
        )
        # compute its square before averaging over realizations
        post_selec["Jz fluctu^2"] = post_selec["Jz fluctu"] ** 2
        res = post_selec.groupby(
            [self.var1.name, self.var2.name], as_index=False
        ).mean()
        # get the count also
        counts = post_selec.groupby(
            [self.var1.name, self.var2.name], as_index=False
        ).count()
        # compute the quantity of interest
        res["denis"] = (res["Jz fluctu^2"] - res["Jz fluctu"] ** 2) / res["N_1*N_2"]
        res["Gab"] = res["N_1*N_2"]  # change name
        res["denis counts"] = counts["Jz fluctu"]

        res_to_keep = res[["Vz1", "Vz2", "Gab", "denis", "denis counts"]]
        self.result = pd.merge(self.result, res_to_keep, on=["Vz1", "Vz2"])
        self.result["P_denis"] = self.result["denis"]
        self.result["log(denis)"] = np.log(self.result["denis"])

    def save_copy_of_total(self):
        """Save a copy of the total dataframe. Important to do if one does bottstraping."""
        self.total_dataframe_copy = copy.deepcopy(self.total)
        self.total_dataframe_copy["Original Cycle"] = self.total_dataframe_copy["Cycle"]
        self.cycles_array_copy = copy.deepcopy(self.cycles_array)
        self.is_there_a_copy_of_total = True

    def recover_true_total(self):
        self.total = copy.deepcopy(self.total_dataframe_copy)
        self.cycles_array = copy.deepcopy(self.cycles_array_copy)

    def bootstrap_total(self):
        """ bootstrap the dataframe total in an efficient way. See 03/01/24 for details.
        We get a matrix from the total dataframe and then we bootstrap it as it is MUCH MORE faster than with pandas.
        """
        if self.is_there_a_copy_of_total is False:
            self.save_copy_of_total()
            print(
                "[Warning] : I just saved a copy of the total dataframe because you will destruct your original dataframe."
            )
        ordata = self.total_dataframe_copy.set_index(
            ["Cycle", self.total_dataframe_copy.groupby("Cycle").cumcount()]
        )
        # ordata.index.names = ["Cycle", "My_tmp_index"]
        original_data_array = ordata.values.reshape(
            (self.n_cycles, -1, len(ordata.columns))
        )
        _, n_tmp_index, _ = original_data_array.shape
        new_indices = np.random.randint(0, self.n_cycles, self.n_cycles)
        NEW = np.zeros_like(original_data_array)
        NEW[:] = original_data_array[new_indices]
        self.total = pd.DataFrame(
            data=NEW.reshape((-1, len(ordata.columns))), columns=list(ordata.columns)
        )
        self.total["Cycle"] = np.repeat(self.cycles_array_copy, n_tmp_index)
        # self.total["My_tmp_index"] = np.tile(np.arange(n_tmp_index), n_cycles)

    ###########################################################
    ############### Plot Functions  #####################
    ###########################################################

    def show_densities(self, cmap="Blues", bins=100, return_fig=False, **kwargs):
        speeds = [("Vx", "Vy"), ("Vx", "Vy"), ("Vx", "Vz"), ("Vy", "Vz")]
        fig, axes = plt.subplots(figsize=(16, 4), ncols=4)

        def draw_box(ax, cX, σX, cY, σY, color="orange", label="box"):
            ax.plot(
                [cX - σX, cX + σX, cX + σX, cX - σX, cX - σX],
                [cY - σY, cY - σY, cY + σY, cY + σY, cY - σY],
                color,
                label=label,
            )

        for i, (nameX, nameY) in enumerate(speeds):
            ax = axes.flatten()[i]
            if i == 0:
                at = self.atoms[self.atoms["Vz"] < 0]
                ax.set_title("Vz < 0")
            elif i == 1:
                at = self.atoms[self.atoms["Vz"] > 0]
                ax.set_title("Vz > 0")
            else:
                at = self.atoms
            X_list = at[nameX].to_numpy()
            Y_list = at[nameY].to_numpy()
            heatmap = ax.hist2d(X_list, Y_list, bins=bins, cmap=cmap, **kwargs)
            for box_num in ["1", "2"]:
                posX = self.boxes[box_num][nameX]["position"]
                sizeX = self.boxes[box_num][nameX]["size"]
                posY = self.boxes[box_num][nameY]["position"]
                sizeY = self.boxes[box_num][nameY]["size"]
                draw_box(
                    ax, posX, sizeX / 2, posY, sizeY / 2, color="orange", label="box"
                )
                draw_box(
                    ax,
                    posX,
                    3 * sizeX / 2,
                    posY,
                    3 * sizeY / 2,
                    color="darkred",
                    label="3*box",
                )
            # ax.legend()
            ax.set_xlabel(nameX)
            ax.set_ylabel(nameY)
            colorbar = plt.colorbar(
                heatmap[3], ax=ax
            )  # Utilisez la quatrième valeur de retour de hist2d

        plt.tight_layout()
        if return_fig:
            return fig
        plt.show()

    def show_density(
        self,
        nameX="Vy",
        nameY="Vz",
        x_bins=100,
        y_bins=100,
        title=None,
        show_boxes=True,
        save=None,
        show_plot=True,
        ax=None,
    ):
        """
            Affiche l'histogramme 2D de colonnes nameX et nameY du dataframe.

        Parameters
        ----------
        nameX, nameY : nom des colonnes du dataframe dont on veut tracer l'histogramme (abscisse, ordonnée) (soit Vx, Vy ou Vz)
        x_bins, y_bins : nombre de bins en abscisse et en ordonnée
        title : titre de l'histogramme
        show_plot : boolean, if we want to see the plot
        save : pathe to save the picture

        Returns
        ----------
        hist_values : The bi-dimensional histogram of samples x and y. Values in x are histogrammed along the first dimension and values in y are histogrammed along the second dimension.
        X_values : The bin edges along the x axis.
        Y_values : The bin edges along the y axis.
        """
        if title == None:
            title = "Histogramme pour {} fichiers".format(self.n_cycles)
        plt.clf()
        # plt.figure(figsize=(10, 7))
        if ax is None:
            fig, ax = plt.subplots()

        X_list = self.atoms[nameX].to_numpy()
        Y_list = self.atoms[nameY].to_numpy()
        hist_values, X_values, Y_values, _ = ax.hist2d(
            X_list, Y_list, bins=[x_bins, y_bins], cmap=plt.cm.Blues
        )
        if show_boxes:

            def draw_box(cX, σX, cY, σY, color="orange", label="box"):
                ax.plot(
                    [cX - σX, cX + σX, cX + σX, cX - σX, cX - σX],
                    [cY - σY, cY - σY, cY + σY, cY + σY, cY - σY],
                    color,
                    label=label,
                )

            # on affiche la boite 1
            posX = self.boxes["1"][nameX]["position"]
            sizeX = self.boxes["1"][nameX]["size"]
            posY = self.boxes["1"][nameY]["position"]
            sizeY = self.boxes["1"][nameY]["size"]
            draw_box(posX, sizeX / 2, posY, sizeY / 2, color="orange", label=None)
            draw_box(
                posX, 3 * sizeX / 2, posY, 3 * sizeY / 2, color="darkred", label=None
            )

            # et le boite 2
            posX = self.boxes["2"][nameX]["position"]
            sizeX = self.boxes["2"][nameX]["size"]
            posY = self.boxes["2"][nameY]["position"]
            sizeY = self.boxes["2"][nameY]["size"]
            draw_box(posX, sizeX / 2, posY, sizeY / 2, color="orange", label="box")
            draw_box(
                posX, 3 * sizeX / 2, posY, 3 * sizeY / 2, color="darkred", label="3*box"
            )
        # ax.colorbar()
        ax.legend()
        ax.set_title(title)
        ax.set_xlabel(nameX)
        ax.set_ylabel(nameY)
        if save:
            fig.savefig(save)
        if show_plot is True:
            plt.show()

        return (hist_values, X_values, Y_values)

    def get_atoms_distribution(self, nbMax, nbPt, posZ, sizeZ, posX, sizeX, posY, sizeY):
        """
        Allows to plot the distribution of the number of atoms averaged over all cycles, either in a box or in an average of boxes
        Parameters
        ----------
        posZ, sizeZ, posX, sizeX, posY, sizeY: parameters of the initial box
        nbMax: number of simultaneous atoms detected to consider
        nbPt: number of boxes on which we average (the 1st box is the one centered on PosZ and of size sizeZ, then the other boxes are the neighbors according to Vz in the direction of increasing Vz and of size sizeZ

        """
        n_cycles = self.n_cycles
        tbin = np.arange(0, nbMax + 1, 1)
        pro = np.zeros(nbMax + 1)
        nb = 0

        for i in range(nbPt):
            box_proba = {
                "Vz": {"position": posZ + sizeZ * i, "size": sizeZ},
                "Vy": {"position": posY, "size": sizeY},
                "Vx": {"position": posX, "size": sizeX},
            }
            atoms_in_box = self.obtain_number_of_atoms_per_cycle_in_box(
                self.atoms, box_proba, "toto"
            )
            nb_bin = int(np.max(atoms_in_box["toto"]))
            a, b = np.histogram(atoms_in_box["toto"], bins=nbMax + 1, range=(0, nbMax))
            pro = a + pro
            if nb_bin > nb:
                nb = nb_bin
        # pro contient la somme des histogrammes sur tous les cycles et toutes les boîtes -> on normalise
        proF = [px / n_cycles / nbPt for px in pro]
        proF = proF[0 : (nb + 1)]
        tbin = tbin[0 : (nb + 1)]
        moy = np.sum(tbin * proF)
        print("population moyenne = ", "{:.3f}".format(moy))
        print("Vérification normalisation : ", "{:.3f}".format(np.sum(proF)))
        yerr = [np.sqrt(px) / n_cycles / nbPt for px in pro]
        yerr = yerr[0 : (nb + 1)]
        # calcul des distributions théoriques
        therm = moy**tbin / (1 + moy) ** (tbin + 1)
        pois = np.exp(-moy) * moy**tbin / factorial(tbin)

        # tracé du graphe
        plt.figure()
        plt.errorbar(tbin, proF, yerr=yerr, label="exp")
        plt.plot(tbin, therm, label="therm")
        plt.plot(tbin, pois, label="pois")
        plt.yscale("log")
        plt.ylim([1e-6, 1])
        plt.xlabel("Number of atoms")
        plt.ylabel("Probability")
        plt.grid(True)
        plt.legend()
        plt.show()

    def plot_population_in_box(self, vz_list_1, sizeZ, posX, sizeX, posY, sizeY):
        """
        Permet de tracer le nombre moyen d'atomes par boîte en fonction du centre de la boîte
        Parameters
        ----------
        vz_list_1 : liste (ou numpy array) des centres à considérer
        sizeZ, posX, sizeX, posY, sizeY : paramètres de la boîte
        """
        nb_atoms_in_box_1 = []

        for k in range(len(vz_list_1)):
            box_proba = {
                "Vz": {"position": vz_list_1[k], "size": sizeZ},
                "Vy": {"position": posY, "size": sizeX},
                "Vx": {"position": posX, "size": sizeY},
            }
            atoms_in_box = self.obtain_number_of_atoms_per_cycle_in_box(
                self.atoms, box_proba, "toto"
            )
            nb_atoms_in_box_1.append(np.mean(atoms_in_box["toto"]))

        plt.figure()
        plt.plot(vz_list_1, nb_atoms_in_box_1)
        plt.grid(True)
        # plt.legend()
        plt.xlabel("center Vz of the box (mm/s)")
        plt.ylabel("mean number of atoms in box")
        plt.show()

    def get_joint_atoms_distribution_in_two_boxes(
        self,
        nbMax,
        nbPt,
        posZ1,
        posZ2,
        sizeZ,
        posX,
        sizeX,
        posY,
        sizeY,
        show=True,
        sizeZ1=None,
        sizeZ2=None,
    ):
        """
        Permet de tracer la distribution jointe moyenne du nombre d'atomes sur tous les cycles, dans deux boîtes données.
        Le nombre de boîtes considérées est impair de façon à toujours prendre autant de boîtes de part et d'autre de la boîte centrale
        ----------
        Parameters
        ----------
        posZ1, posZ2 : centres des boîtes à scanner
        sizeZ, posX, sizeX, posY, sizeY : paramètres de la boîte à scanner
        nbMax : nombre d'atomes simultanés détectés à considérer
        nbPt : nombre de boîtes sur lesquelles on moyenne (la 1e boîte est celle centrée sur PosZ et de taille sizeZ, puis les autres boîtes sont les voisines selon Vz dans le sens des Vz croissants et de taille sizeZ
        /!\ nbPt doit être impair
        show = est ce que on montre la 2D map
        """
        n_cycles = self.n_cycles
        if sizeZ1 == sizeZ2:
            if type(sizeZ1) == type(None):
                sizeZ1 = sizeZ
                sizeZ2 = sizeZ

        # fonction auxiliaire pour construire les boîtes
        def return_boxes(
            i,
            posZ1,
            posZ2,
            posX,
            posY,
            sizeZ,
            sizeX,
            sizeY,
            sizeZ1=sizeZ1,
            sizeZ2=sizeZ2,
        ):

            assert nbPt % 2 == 1
            if i == 0:
                box1_proba = {
                    "Vz": {"position": posZ1, "size": sizeZ1},
                    "Vy": {"position": posY, "size": sizeY},
                    "Vx": {"position": posX, "size": sizeX},
                }
                box2_proba = {
                    "Vz": {"position": posZ2, "size": sizeZ1},
                    "Vy": {"position": posY, "size": sizeY},
                    "Vx": {"position": posX, "size": sizeX},
                }
            else:
                if i % 2 == 1:
                    box1_proba = {
                        "Vz": {"position": posZ1 + sizeZ * i, "size": sizeZ},
                        "Vy": {"position": posY, "size": sizeY},
                        "Vx": {"position": posX, "size": sizeX},
                    }
                    box2_proba = {
                        "Vz": {"position": posZ2 - sizeZ * i, "size": sizeZ},
                        "Vy": {"position": posY, "size": sizeY},
                        "Vx": {"position": posX, "size": sizeX},
                    }
                else:
                    box1_proba = {
                        "Vz": {"position": posZ1 - sizeZ * i, "size": sizeZ},
                        "Vy": {"position": posY, "size": sizeY},
                        "Vx": {"position": posX, "size": sizeX},
                    }
                    box2_proba = {
                        "Vz": {"position": posZ2 + sizeZ * i, "size": sizeZ},
                        "Vy": {"position": posY, "size": sizeY},
                        "Vx": {"position": posX, "size": sizeX},
                    }
            return (box1_proba, box2_proba)

        tbin = np.arange(0, nbMax + 1, 1)
        pro2D = np.zeros([nbMax + 1, nbMax + 1], dtype=int)
        pro1 = np.zeros(nbMax + 1, dtype=int)
        pro2 = np.zeros(nbMax + 1, dtype=int)
        nb = 0
        for i in range(nbPt):
            box1_proba, box2_proba = return_boxes(
                i, posZ1, posZ2, posX, posY, sizeZ, sizeX, sizeY
            )
            atoms1_in_box = self.obtain_number_of_atoms_per_cycle_in_box(
                self.atoms, box1_proba, "toto"
            )
            atoms2_in_box = self.obtain_number_of_atoms_per_cycle_in_box(
                self.atoms, box2_proba, "toto"
            )
            a2D, b, c = np.histogram2d(
                atoms1_in_box["toto"],
                atoms2_in_box["toto"],
                bins=nbMax + 1,
                range=([[0, nbMax], [0, nbMax]]),
            )
            a1, b = np.histogram(
                atoms1_in_box["toto"], bins=nbMax + 1, range=(0, nbMax)
            )
            a2, b = np.histogram(
                atoms2_in_box["toto"], bins=nbMax + 1, range=(0, nbMax)
            )
            pro2D = a2D + pro2D
            pro1 = a1 + pro1
            pro2 = a2 + pro2
        moy1 = np.sum(tbin * pro1) / n_cycles / nbPt
        moy2 = np.sum(tbin * pro2) / n_cycles / nbPt
        print(
            "Population moyenne zone 1 : ",
            "{:.3f}".format(moy1),
            "  ---   Population moyenne zone 2 : ",
            "{:.3f}".format(moy2),
        )
        pro2D = pro2D / nbPt / n_cycles
        if show:
            im = plt.imshow(
                pro2D,
                interpolation="nearest",
                origin="lower",
                cmap=cm.rainbow,
                norm=colors.LogNorm(),
            )
            plt.colorbar(im)
            plt.show()
        return (moy1, moy2, pro2D)


class Variable:
    """
    Defines the Variable class. This class contains the information of the box to count/integrate atoms and compute correlations
    """

    def __init__(self, **kwargs):
        self.box = "1"  # the box number of the scanned parameter (1 or 2)
        self.axe = "Vx"  # the axis of the scanned parameter (Vx, Vy or Vz)
        self.type = "size"  # type: size or position, that is if we scan size of box or position
        self.name = "ΔVx"  # its name for the column in the dataframe

        """ The scan can be defined in two ways: either we give a list of values or we provide min, max and step for scan"""
        # min, max and step of scan
        self.min = 0  # start of scan
        self.max = 2 # end of scan
        self.step = 4 # set for scan
        # list of points for scan        
        self.values = []  # list of points
        self.round_decimal = 7
        self.__dict__.update(kwargs)
        # if no list is provided, we create the list based on min, max and step defined by user
        if len(self.values) == 0:
            self.built_values() # create list
        # get new min, nax and step in case they are not defined or numpy does something funny
        self.get_values_caracteristics()

    def built_values(self):
        """ Creates list of values to scan using min, max and step defined by user """
        mini = min(self.min, self.max)
        maxi = max(self.min, self.max)
        self.values = np.arange(mini, maxi, self.step)
        self.values = np.round(self.values, self.round_decimal)

    def built_name(self, add_box_number=True):
        """built a default name"""
        name = ""
        if self.type == "size":
            name += "Δ"
        name += self.axe
        if add_box_number:
            name += self.box
        return name

    def get_values_caracteristics(self):
        """ gets min, max and step of self.values """
        self.values = np.array(self.values)
        self.min = np.min(self.values)
        self.max = np.max(self.values)
        self.n_step = len(self.values)

    def get_value_i(self, i):
        """
        Returns the i-th value of self.values
        """
        return self.values[i]


class CorrelationXYIntegrated(Correlation):
    def __init__(self, atoms, NsliceX, NsliceY, **kwargs):
        self.NsliceX = self.check_Nslices_value(NsliceX)
        self.NsliceY = self.check_Nslices_value(NsliceY)
        super().__init__(atoms, **kwargs)
        self.global_boxes = copy.deepcopy(self.boxes)
        self.built_boxes_list()

    def check_Nslices_value(self, N):
        if N < 1:
            print(
                "WARNING : Nx or Ny value seems strange in CorrelationXYIntegrated Initialisation. Setting to 1."
            )
            return 1
        if int(N) != N:
            print("WARNING : Nx or Ny value does not seem an integer. Setting to 1.")
            return round(N)
        return N

    def compute_correlations_XYintegrated(self):
        all_result = []
        for idx, box in enumerate(self.boxes_list):
            self.set_boxes(box)
            self.compute_correlations()
            all_result.append(self.result)
        self.all_result = pd.concat(all_result)
        self.all_result["<N1><N2>"] = self.all_result["N_1"] * self.all_result["N_2"]
        self.integrated_result = self.all_result.groupby(
            [self.var1.name, self.var2.name], as_index=False
        ).sum()
        self.integrated_result["g^2"] = (
            self.integrated_result["N_1*N_2"] / self.integrated_result["<N1><N2>"]
        )
        self.integrated_result["normalized variance"] = self.integrated_result[
            "variance"
        ] / (self.integrated_result["N_1"] + self.integrated_result["N_2"])
        # je remets la boite initial.
        self.set_boxes(self.global_boxes)
        self.result = self.integrated_result

    def built_boxes_list(self):
        self.boxes_list = []
        import itertools as itt

        for i, XY in enumerate(
            list(itt.product(np.arange(self.NsliceX), np.arange(self.NsliceY)))
        ):
            self.boxes_list.append(
                {
                    "1": {
                        "Vx": {
                            "size": self.global_boxes["1"]["Vx"]["size"] / self.NsliceX,
                            "position": self.global_boxes["1"]["Vx"]["position"]
                            - self.global_boxes["1"]["Vx"]["size"] / 2.0
                            + self.global_boxes["1"]["Vx"]["size"]
                            * (0.5 + XY[0])
                            / self.NsliceX,
                        },
                        "Vy": {
                            "size": self.global_boxes["1"]["Vy"]["size"] / self.NsliceY,
                            "position": self.global_boxes["1"]["Vy"]["position"]
                            - self.global_boxes["1"]["Vy"]["size"] / 2.0
                            + self.global_boxes["1"]["Vy"]["size"]
                            * (0.5 + XY[1])
                            / self.NsliceY,
                        },
                        "Vz": {
                            "size": self.global_boxes["1"]["Vz"]["size"],
                            "position": self.global_boxes["1"]["Vz"]["position"],
                        },
                    },
                    "2": {
                        "Vx": {
                            "size": self.global_boxes["2"]["Vx"]["size"] / self.NsliceX,
                            "position": self.global_boxes["2"]["Vx"]["position"]
                            - self.global_boxes["2"]["Vx"]["size"] / 2.0
                            + self.global_boxes["2"]["Vx"]["size"]
                            * (0.5 + XY[0])
                            / self.NsliceX,
                        },
                        "Vy": {
                            "size": self.global_boxes["2"]["Vy"]["size"] / self.NsliceY,
                            "position": self.global_boxes["2"]["Vy"]["position"]
                            - self.global_boxes["2"]["Vy"]["size"] / 2.0
                            + self.global_boxes["2"]["Vy"]["size"]
                            * (0.5 + XY[1])
                            / self.NsliceY,
                        },
                        "Vz": {
                            "size": self.global_boxes["2"]["Vz"]["size"],
                            "position": self.global_boxes["2"]["Vz"]["position"],
                        },
                    },
                }
            )


class CorrelationCollision(Correlation):
    def __init__(self, atoms, NsliceTheta, NslicePhi, **kwargs):
        self.NslicePhi = self.check_Nslices_value(NslicePhi)
        self.NsliceTheta = self.check_Nslices_value(NsliceTheta)
        super().__init__(atoms, **kwargs)
        self.compute_spherical_coordinates()
        self.built_boxes_list()

    def check_Nslices_value(self, N):
        if N < 2:
            print(
                "WARNING : NTheta or NPhi value seems strange in CorrelationCollision Initialisation. Setting to 2."
            )
            return 2
        if int(N) != N:
            val = max(round(N), 2)
            print(
                "WARNING :  NTheta or NPhi value does not seem an integer. Setting to {}.".format(
                    val
                )
            )

            return val
        return N

    def built_boxes_list(self):
        theta_min = np.min(self.atoms["theta"])
        theta_max = np.max(self.atoms["theta"])
        phi_min = np.min(self.atoms["phi"])
        phi_max = np.max(self.atoms["phi"])
        liste_theta = np.linspace()


class CorrelationThirdOrder(Correlation):
    var1 = None
    var2 = None
    var3 = None

    def __init__(self, atoms, **kwargs):
        super().__init__(atoms, **kwargs)
        self.boxes = {
            "1": {
                "Vx": {"size": 80, "position": 0},
                "Vy": {"size": 80, "position": 0},
                "Vz": {"size": 0.9, "position": 56},
            },
            "2": {
                "Vx": {"size": 80, "position": 0},
                "Vy": {"size": 80, "position": 0},
                "Vz": {"size": 0.9, "position": 130},
            },
            "3": {
                "Vx": {"size": 80, "position": 0},
                "Vy": {"size": 80, "position": 0},
                "Vz": {"size": 0.9, "position": 130},
            },
        }
        self.__dict__.update(kwargs)

    def check_variables(self):

        for i, var in enumerate([self.var1, self.var2, self.var3]):
            if var == None:
                log.warning(
                    f"[CorrelationsThirdOrder] Variable {i+1} was not defined. Please define it."
                )
                return 1
        return 0

    def compute_correlations(self):
        if self.check_variables():
            return
        total = self.generate_total()
        self.compute_result(total)

    def define_variable3(self, **kwargs):
        self.var3 = Variable(**kwargs)

    def generate_total(self):
        """method that generates the total dataframe to compute correlation.
        It gets the the number of atoms per position/size of the 3 scanned variables. Then, it computes the dataframe which is the cartesian product of the three.
        """
        ## the following is just a copypaste of the compute_correlations_different_box_scanned method from Correlation
        box = self.boxes[self.var1.box].copy()
        posi_and_size = box.pop(self.var1.axe)
        scanned_box = {self.var1.axe: posi_and_size}
        df_atoms_var1 = self.get_atoms_in_box(self.atoms, box)
        # result_var1 is a table with the cycle, the number of atoms in the cycle and the variable 1 value.
        result_var1 = self.counts_atoms_in_boxes_one_variable(
            df_atoms_var1, self.var1, scanned_box, column_name="N_" + self.var1.box
        )
        # -- Do the same for the second variable
        box = self.boxes[self.var2.box].copy()
        posi_and_size = box.pop(self.var2.axe)
        scanned_box = {self.var2.axe: posi_and_size}
        df_atoms_var2 = self.get_atoms_in_box(self.atoms, box)
        result_var2 = self.counts_atoms_in_boxes_one_variable(
            df_atoms_var2, self.var2, scanned_box, column_name="N_" + self.var2.box
        )
        # -- idem with the third variable
        box = self.boxes[self.var3.box].copy()
        posi_and_size = box.pop(self.var3.axe)
        scanned_box = {self.var3.axe: posi_and_size}
        df_atoms_var3 = self.get_atoms_in_box(self.atoms, box)
        result_var3 = self.counts_atoms_in_boxes_one_variable(
            df_atoms_var3, self.var3, scanned_box, column_name="N_" + self.var3.box
        )

        # %#% STEP2
        # On construit le dataframe total, qui initialement contient 5 colonnes : Cycle le cycle, N_1 et N_2 nombre d'atomes dans la boîte 1 et 2, self.var1 et self.var2 la position/taille des boîtes lors du scan. Le nombre de lignes de total est dont Nombre_de_cycles x Nombre_de_différentes_var1 x Nombre_de_différentes_var2.
        total = pd.merge(result_var1, result_var2, on="Cycle")
        total = pd.merge(total, result_var3)
        return total

    def compute_result(self, total):
        total["N1N2"] = total["N_1"] * total["N_2"]
        total["N1N3"] = total["N_1"] * total["N_3"]
        total["N2N3"] = total["N_2"] * total["N_3"]
        total["N1N2N3"] = total["N_1"] * total["N_2"] * total["N_3"]
        self.total = total
        self.result = total.groupby(
            [self.var1.name, self.var2.name, self.var3.name], as_index=False
        ).mean()
        ## - define connected correlations
        self.result["Gc12"] = (
            self.result["N1N2"] - self.result["N_1"] * self.result["N_2"]
        )
        self.result["Gc13"] = (
            self.result["N1N3"] - self.result["N_1"] * self.result["N_3"]
        )
        self.result["Gc23"] = (
            self.result["N2N3"] - self.result["N_2"] * self.result["N_3"]
        )
        self.result["Gc123"] = (
            self.result["N1N2N3"]
            - (
                self.result["N_3"] * self.result["Gc12"]
                + self.result["N_2"] * self.result["Gc13"]
                + self.result["N_1"] * self.result["Gc23"]
            )
            - self.result["N_1"] * self.result["N_2"] * self.result["N_3"]
        )

        self.result["g2_12"] = self.result["N1N2"] / (
            self.result["N_1"] * self.result["N_2"]
        )
        self.result["g2_23"] = self.result["N2N3"] / (
            self.result["N_3"] * self.result["N_2"]
        )
        self.result["g2_13"] = self.result["N1N3"] / (
            self.result["N_1"] * self.result["N_3"]
        )


class CorrelationFourthOrder(Correlation):
    var1 = None
    var2 = None
    var3 = None
    var4 = None

    def __init__(self, atoms, **kwargs):
        super().__init__(atoms, **kwargs)
        self.boxes = {
            "1": {
                "Vx": {"size": 80, "position": 0},
                "Vy": {"size": 80, "position": 0},
                "Vz": {"size": 0.9, "position": 56},
            },
            "2": {
                "Vx": {"size": 80, "position": 0},
                "Vy": {"size": 80, "position": 0},
                "Vz": {"size": 0.9, "position": 130},
            },
            "3": {
                "Vx": {"size": 80, "position": 0},
                "Vy": {"size": 80, "position": 0},
                "Vz": {"size": 0.9, "position": 130},
            },
            "4": {
                "Vx": {"size": 80, "position": 0},
                "Vy": {"size": 80, "position": 0},
                "Vz": {"size": 0.9, "position": 130},
            },
        }
        self.__dict__.update(kwargs)

    def check_variables(self):

        for i, var in enumerate([self.var1, self.var2, self.var3, self.var4]):
            if var == None:
                log.warning(
                    f"[CorrelationsThirdOrder] Variable {i+1} was not defined. Please define it."
                )
                return 1
        return 0

    def compute_correlations(self):
        if self.check_variables():
            return
        total = self.generate_total()
        self.compute_result(total)

    def define_variable3(self, **kwargs):
        self.var3 = Variable(**kwargs)

    def define_variable4(self, **kwargs):
        self.var4 = Variable(**kwargs)

    def generate_total(self):
        """method that generates the total dataframe to compute correlation.
        It gets the the number of atoms per position/size of the 3 scanned variables. Then, it computes the dataframe which is the cartesian product of the three.
        """
        ## the following is just a copypaste of the compute_correlations_different_box_scanned method from Correlation
        box = self.boxes[self.var1.box].copy()
        posi_and_size = box.pop(self.var1.axe)
        scanned_box = {self.var1.axe: posi_and_size}
        df_atoms_var1 = self.get_atoms_in_box(self.atoms, box)
        # result_var1 is a table with the cycle, the number of atoms in the cycle and the variable 1 value.
        result_var1 = self.counts_atoms_in_boxes_one_variable(
            df_atoms_var1, self.var1, scanned_box, column_name="N_" + self.var1.box
        )
        # -- Do the same for the second variable
        box = self.boxes[self.var2.box].copy()
        posi_and_size = box.pop(self.var2.axe)
        scanned_box = {self.var2.axe: posi_and_size}
        df_atoms_var2 = self.get_atoms_in_box(self.atoms, box)
        result_var2 = self.counts_atoms_in_boxes_one_variable(
            df_atoms_var2, self.var2, scanned_box, column_name="N_" + self.var2.box
        )
        at_random = result_var2["N_" + self.var2.box].to_numpy()
        np.random.shuffle(at_random)
        # result_var2["N_" + self.var2.box] = at_random
        # -- idem with the third variable
        box = self.boxes[self.var3.box].copy()
        posi_and_size = box.pop(self.var3.axe)
        scanned_box = {self.var3.axe: posi_and_size}
        df_atoms_var3 = self.get_atoms_in_box(self.atoms, box)
        result_var3 = self.counts_atoms_in_boxes_one_variable(
            df_atoms_var3, self.var3, scanned_box, column_name="N_" + self.var3.box
        )
        at_random = result_var3["N_" + self.var3.box].to_numpy()
        np.random.shuffle(at_random)
        # result_var3["N_" + self.var3.box] = at_random

        # -- idem with the fourth variable
        box = self.boxes[self.var4.box].copy()
        posi_and_size = box.pop(self.var4.axe)
        scanned_box = {self.var4.axe: posi_and_size}
        df_atoms_var4 = self.get_atoms_in_box(self.atoms, box)
        result_var4 = self.counts_atoms_in_boxes_one_variable(
            df_atoms_var4, self.var4, scanned_box, column_name="N_" + self.var4.box
        )
        at_random = result_var4["N_" + self.var4.box].to_numpy()
        np.random.shuffle(at_random)
        # result_var4["N_" + self.var4.box] = at_random

        # %#% STEP2
        # On construit le dataframe total, qui initialement contient 5 colonnes : Cycle le cycle, N_1 et N_2 nombre d'atomes dans la boîte 1 et 2, self.var1 et self.var2 la position/taille des boîtes lors du scan. Le nombre de lignes de total est dont Nombre_de_cycles x Nombre_de_différentes_var1 x Nombre_de_différentes_var2.
        total = pd.merge(result_var1, result_var2, on="Cycle")
        total = pd.merge(total, result_var3, on="Cycle")
        total = pd.merge(total, result_var4, on="Cycle")
        return total

    def compute_result(self, total):
        total["N1N2"] = total["N_1"] * total["N_2"]
        total["N1N3"] = total["N_1"] * total["N_3"]
        total["N1N4"] = total["N_1"] * total["N_4"]
        total["N2N3"] = total["N_2"] * total["N_3"]
        total["N2N4"] = total["N_2"] * total["N_4"]
        total["N3N4"] = total["N_3"] * total["N_4"]
        total["N1N2N3"] = total["N_1"] * total["N_2"] * total["N_3"]
        total["N1N2N4"] = total["N_1"] * total["N_2"] * total["N_4"]
        total["N1N3N4"] = total["N_1"] * total["N_3"] * total["N_4"]
        total["N2N3N4"] = total["N_2"] * total["N_3"] * total["N_4"]
        total["N1N2N3N4"] = total["N_1"] * total["N_2"] * total["N_3"] * total["N_4"]
        self.total = total
        self.result = total.groupby(
            [self.var1.name, self.var2.name, self.var3.name, self.var4.name],
            as_index=False,
        ).mean()
        ## - define connected correlations
        self.result["Gc12"] = (
            self.result["N1N2"] - self.result["N_1"] * self.result["N_2"]
        )
        self.result["Gc13"] = (
            self.result["N1N3"] - self.result["N_1"] * self.result["N_3"]
        )
        self.result["Gc14"] = (
            self.result["N1N4"] - self.result["N_1"] * self.result["N_4"]
        )
        self.result["Gc23"] = (
            self.result["N2N3"] - self.result["N_2"] * self.result["N_3"]
        )
        self.result["Gc24"] = (
            self.result["N2N4"] - self.result["N_2"] * self.result["N_4"]
        )
        self.result["Gc34"] = (
            self.result["N3N4"] - self.result["N_3"] * self.result["N_4"]
        )
        self.result["Gc123"] = (
            self.result["N1N2N3"]
            - self.result["Gc12"]
            - self.result["Gc13"]
            - self.result["Gc23"]
        )
        self.result["Gc124"] = (
            self.result["N1N2N4"]
            - self.result["Gc12"]
            - self.result["Gc14"]
            - self.result["Gc24"]
        )
        self.result["Gc134"] = (
            self.result["N1N3N4"]
            - self.result["Gc13"]
            - self.result["Gc14"]
            - self.result["Gc34"]
        )
        self.result["Gc234"] = (
            self.result["N2N3N4"]
            - self.result["Gc23"]
            - self.result["Gc24"]
            - self.result["Gc34"]
        )
        ## - define connected correlations
        self.result["Gc12"] = (
            self.result["N1N2"] - self.result["N_1"] * self.result["N_2"]
        )
        self.result["Gc13"] = (
            self.result["N1N3"] - self.result["N_1"] * self.result["N_3"]
        )
        self.result["Gc23"] = (
            self.result["N2N3"] - self.result["N_2"] * self.result["N_3"]
        )
        self.result["Gc123"] = (
            self.result["N1N2N3"]
            - (
                self.result["N_3"] * self.result["Gc12"]
                + self.result["N_2"] * self.result["Gc13"]
                + self.result["N_1"] * self.result["Gc23"]
            )
            - self.result["N_1"] * self.result["N_2"] * self.result["N_3"]
        )
        ### - define fourth order correlations
        self.result["Gc1234"] = (
            self.result["N1N2N3N4"]
            - (
                self.result["N_1"] * self.result["N2N3N4"]
                + self.result["N_2"] * self.result["N1N3N4"]
                + self.result["N_3"] * self.result["N1N2N4"]
                + self.result["N_4"] * self.result["N1N2N3"]
                + self.result["N1N2"] * self.result["N3N4"]
                + self.result["N1N3"] * self.result["N2N4"]
                + self.result["N1N4"] * self.result["N2N3"]
            )
            + 2
            * (
                self.result["N1N2"] * self.result["N_3"] * self.result["N_4"]
                + self.result["N1N3"] * self.result["N_2"] * self.result["N_4"]
                + self.result["N1N4"] * self.result["N_2"] * self.result["N_3"]
                + self.result["N2N3"] * self.result["N_1"] * self.result["N_4"]
            )
        )
        -6 * self.result["N_1"] * self.result["N_2"] * self.result["N_3"] * self.result[
            "N_4"
        ]
        self.result["gc1234"] = self.result["Gc1234"] / (
            self.result["N_1"]
            * self.result["N_2"]
            * self.result["N_3"]
            * self.result["N_4"]
        )
        self.result["g2_12"] = self.result["N1N2"] / (
            self.result["N_1"] * self.result["N_2"]
        )
        self.result["g2_23"] = self.result["N2N3"] / (
            self.result["N_3"] * self.result["N_2"]
        )
        self.result["g2_13"] = self.result["N1N3"] / (
            self.result["N_1"] * self.result["N_3"]
        )


if __name__ == "__main__":
    import seaborn as sns
    import pandas as pd

    selected_data = pd.read_pickle(
        "/home/victor/ownCloud/LabWiki/Journal/2023/06/13/data.pkl"
    )
    selec_bec_arrival_times = pd.read_pickle(
        "/home/victor/ownCloud/LabWiki/Journal/2023/06/13/bec.pkl"
    )
    ROI = {
        "Vz": {"min": -50, "max": 50},
        "Vy": {"min": -70, "max": 70},
        "Vx": {"max": 25, "min": -75},
    }
    boxZsize = 5
    boxXsize = 10
    boxYsize = 10
    Xposition = 0
    Yposition = 0
    boxes = {
        "1": {
            "Vx": {"size": boxXsize, "position": Xposition},
            "Vy": {"size": boxYsize, "position": Yposition},
            "Vz": {"size": boxZsize, "position": 25},
        },
        "2": {
            "Vx": {"size": boxXsize, "position": Xposition},
            "Vy": {"size": boxYsize, "position": Yposition},
            "Vz": {"size": boxZsize, "position": -25},
        },
    }
    corr = Correlation(
        selected_data,
        ROI=ROI,
        boxes=boxes,
        raman_kick=42.5,
        bec_arrival_time=selec_bec_arrival_times["BEC Arrival Time"].mean(),
        ref_frame_speed={"Vx": -2, "Vy": -5, "Vz": 94},
        remove_shot_noise=False,
    )
    corr.define_variable1(
        box="1", axe="Vx", type="position", name="Vx1", min=-20, max=11, step=10
    )
    corr.define_variable2(
        box="1", axe="Vy", type="position", name="Vy1", min=-20, max=20, step=10
    )
    corr.compute_correlations()
    df_pivoted_correlations = corr.result.pivot(
        index="Vx1", columns="Vy1", values="g^2"
    )
    fig, axes = plt.subplots(figsize=(8, 3), ncols=2)
    sns.heatmap(
        df_pivoted_correlations,
        cmap="seismic",
        ax=axes[0],
        # norm=LogNorm()
        center=1,
    )
    axes[0].invert_yaxis()
    sns.scatterplot(data=corr.result, x="Vx1", y="N_1", hue="Vy1", palette="Dark2")
    plt.tight_layout()
    plt.show()
