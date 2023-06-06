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

Definition of Correlation and Variable classes

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


class Correlation:
    """
    Classe corrélation. Prend en entrée un dataframe avec X, Y et T et calcule de corrélations etc...

    Mandatory parameters
    --------------------
    atoms : pandas dataframe de 4 colonnes comportant le numéro de cycle, et les temps et positions d'arrivée des atomes. Les colonnes devront s'appeler "Cycle" ; "X" ; "Y" et "T".

    n_cycles : le nombre de fichiers / cycles.

    Main attributs
    --------------------
    atoms : pandas dataframe de 4 colonnes comportant le numéro de cycle, et les vitesses des atomes selon les 3 axes Vx, Vy et Vz



    Other attributs
    --------------------
    boxes : dictionnaire donnant position et taille des boîtes 1 et 2 pour les corrélations.
        Forme : { "1": {"Vx": {"size": 10, "position": 0},
                        "Vy": {"size": 3, "position": 0},
                        "Vz": {"size": 0.9, "position": 56},   },
                "2": {  "Vx": {"size": 10, "position": 0},
                        "Vy": {"size": 3, "position": 0},
                        "Vz": {"size": 0.9, "position": 80},   },}

    ROI : dictionnaire, Region Of Interest
        définir une ROI adaptée permet de réduire la taille des dataframe utilisés (en mm/s)
        ATTENTION : l'application de la méthode apply_ROI() agit sur le dataframe atoms donc si vous vous rendez compte que votre ROI est trop petite, il faut recharger le dataframe atoms.
        Format : {"Vx": {"max":120, "min":-120}}

    ROD : dictionnaire, Region Of Desinterest
        pour ne sélectionner que les atomes en dehors de la ROD. Même remarque que pour ROI et même format.

    bec_arrival_time : float, temps d'arrivée du condensat en ms

    raman_kick : float, kick raman en mm/s

    var1 et var2 : objet Variable (cf classe ci-dessous), les paramètres des boîtes que nous allons changer pour faire les corrélations.

    round_decimal : il s'est avéré (LabJournal du 24/05/2022) que python fait des arrondis un peu bizarre lorsqu'il calcule Vz1 + Vz2 : j'arrondi donc tous les nombres concernant Vz1 et Vz2 (ou plutot les varaibels self.var1.name) à la décimale  round_decimal ( par défaut 5) (--> voir la méthode copute_result)


    Some Methods
    --------------------
    apply_ROI() / apply_ROD() : select atoms only in ROI / outside ROD.

    define_variable1 / define_variable2 : constructs the variable on which we will compute corrrelations. See Variable class for more infos.

    compute_correlations : construct the result dataframe in which correlations are stored.

    show_density : plot the density of the atoms dataframe

    """

    def __init__(self, atoms, **kwargs):
        """
        Object initialization, sets parameters as the user defined, build the atoms dataframe and apply ROD and ROI.
        """
        self.atoms = copy.deepcopy(atoms)
        self.cycles_array = atoms["Cycle"].unique()
        self.n_cycles = len(atoms["Cycle"].unique())
        self.bec_arrival_time = 307.763  # temps d'arrivée du BEC, en ms
        self.theoretical_arrival_time = 307.763  # 24/06/2022 & 17/05/2022
        self.raman_kick = 42.5  # mm/s, kick Raman
        self.gravity = 9.81
        self.bad_shot_limit = 100
        self.var1 = None
        self.var2 = None
        self.remove_shot_noise = True
        self.ref_frame_speed = {"Vx": 0, "Vy": 0, "Vz": 0}
        self.ROI = {}  # Region Of Interest
        self.ROD = {}  # Region Of Desinterest.
        self.round_decimal = 7
        self.id = int(time.time())
        self.is_there_a_copy_of_atoms = False
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
        self.compute_errors = False
        self.__dict__.update(kwargs)
        self.bec_arrival_time = copy.deepcopy(self.bec_arrival_time)
        self.boxes = copy.deepcopy(self.boxes)
        # atoms : (mm,mm, ms) --> (mm/s, mm/s, mm/s)
        self.build_the_atoms_dataframe()
        self.apply_ROI()  # Keep only atoms in ROI
        self.apply_ROD()  # Take off atoms in Region Of Desinterest.

        # print("Data are loaded")
        # Initialisation de self.result datframe avec toutes les corrélations
        # self.result = pd.DataFrame(
        #     np.zeros(1, len(self.quantity_of_interest())),
        #     columns=self.quantity_of_interest,
        # )

    def set_boxes(self, boxes):
        self.boxes = boxes.copy()

    def build_the_atoms_dataframe(self):
        """
        Cette méthode construit le dataframe contenant l'ensemble des positions des atomes : elle construit le dataframe self.atoms, dont les vitesses sont exprimées en mm/s à partir du dataframe initial.
        """
        # Cleaning the self.atoms dataframe : must contains only 4 columns (to avoid mixing when merging.)
        for column in self.atoms.columns:
            if column not in ["X", "Y", "T", "Cycle"]:
                self.atoms.drop(column, inplace=True, axis=1)
        if (
            type(self.bec_arrival_time) == int
            or type(self.bec_arrival_time) == float
            or isinstance(self.bec_arrival_time, np.float64)
            or isinstance(self.bec_arrival_time, np.float32)
            or isinstance(self.bec_arrival_time, np.int32)
            or isinstance(self.bec_arrival_time, np.int64)
        ):
            l_fall = 0.5 * self.gravity * self.bec_arrival_time**2
            self.atoms["X"] = 1000 * self.atoms["X"] / self.atoms["T"] + self.raman_kick
            self.atoms = self.atoms.rename(columns={"X": "Vx"})
            self.atoms["Y"] = 1000 * self.atoms["Y"] / self.atoms["T"]
            self.atoms = self.atoms.rename(columns={"Y": "Vy"})
            self.atoms["T"] = (
                0.5 * self.gravity * self.atoms["T"] - l_fall / self.atoms["T"]
            )
            self.atoms = self.atoms.rename(columns={"T": "Vz"})
        elif isinstance(self.bec_arrival_time, pd.DataFrame):
            self.atoms["X"] = (
                1000 * self.atoms["X"] / self.atoms["T"]
            )  # + self.raman_kick
            self.atoms = self.atoms.rename(columns={"X": "Vx"})
            self.atoms["Y"] = 1000 * self.atoms["Y"] / self.atoms["T"]
            self.atoms = self.atoms.rename(columns={"Y": "Vy"})

            l_fall = 0.5 * self.gravity * self.theoretical_arrival_time**2
            self.atoms["T"] = (
                0.5 * self.gravity * self.atoms["T"] - l_fall / self.atoms["T"]
            )
            self.atoms = self.atoms.rename(columns={"T": "Vz"})
            ## Start to recenter data
            # compute the speed of BECs in mm/s, starting with Vx, Vy and finally Vz.

            if "BEC Center X" in self.bec_arrival_time.columns:
                self.bec_arrival_time["BEC X"] = (
                    1000
                    * self.bec_arrival_time["BEC Center X"]
                    / self.bec_arrival_time["BEC Arrival Time"]
                )
            else:
                self.bec_arrival_time["BEC X"] = -self.raman_kick
            if "BEC Center Y" in self.bec_arrival_time.columns:
                self.bec_arrival_time["BEC Y"] = (
                    1000
                    * self.bec_arrival_time["BEC Center Y"]
                    / self.bec_arrival_time["BEC Arrival Time"]
                )
            else:
                self.bec_arrival_time["BEC Y"] = 0
            self.bec_arrival_time["BEC Z"] = (
                0.5 * self.gravity * self.bec_arrival_time["BEC Arrival Time"]
                - l_fall / self.bec_arrival_time["BEC Arrival Time"]
            )
            # check now if we have perform some Bragg diffraction to fit the BEC
            self.atoms = self.merge_dataframe_on_cycles(
                self.atoms, self.bec_arrival_time
            )
            # take off the speed of each BEC
            self.atoms["Vz"] = self.atoms["Vz"] - self.atoms["BEC Z"]
            self.atoms["Vx"] = self.atoms["Vx"] - self.atoms["BEC X"]
            self.atoms["Vy"] = self.atoms["Vy"] - self.atoms["BEC Y"]
            # print(self.bec_arrival_time["BEC Arrival Time"].mean())

            ## On supprime ensuite les bad shots : là où il y a moins de 100 atomes
            df = self.atoms[self.atoms["Number of Atoms"] < self.bad_shot_limit]
            for cycle, nb_atoms in zip(df["Cycle"], df["Number of Atoms"]):
                print(f"Delete cycle {cycle} with only {nb_atoms} atoms.")
                self.n_cycles -= 1
            # drop columns with no more interest (for later calculations)
            for column in self.atoms.columns:
                if column not in ["Vz", "Vx", "Vy", "Cycle", "Vperp", "theta"]:
                    self.atoms.drop(column, inplace=True, axis=1)

        else:
            print(
                "[ERROR] From build_the_atoms_dataframe : the bec_arrival_time instance is not recognized."
            )

        for axis in self.ref_frame_speed:
            if axis in ["Vx", "Vy", "Vz"] and self.ref_frame_speed[axis] != 0:
                print(
                    f"[INFO] : Reference frame is moving at {self.ref_frame_speed[axis]} mm/s along the {axis} axis."
                )
                self.atoms[axis] -= self.ref_frame_speed[axis]

        self.cycles_array = self.atoms["Cycle"].unique()
        self.n_cycles = len(self.atoms["Cycle"].unique())
        self.compute_cylindrical_coordinates()

    def compute_cylindrical_coordinates(self):
        """Compute transverse velocity and angular angle of the atom dataframe."""
        self.atoms["Vperp"] = np.sqrt(self.atoms["Vx"] ** 2 + self.atoms["Vy"] ** 2)
        self.atoms["theta"] = np.arccos(
            self.atoms["Vx"] / self.atoms["Vperp"]
        )  # [0, pi] range
        local_condition = self.atoms["Vy"] < 0
        self.atoms.loc[local_condition, "theta"] = 2 * np.pi - self.atoms["theta"]

    def save_copy_of_atoms(self):
        """Save a copy of the atom dataframe. Important to do if one does bottstraping."""
        self.atoms_dataframe_copy = copy.deepcopy(self.atoms)
        self.cycles_array_copy = copy.deepcopy(self.cycles_array)
        self.is_there_a_copy_of_atoms = True

    def recover_true_atoms(self):
        self.atoms = copy.deepcopy(self.atoms_dataframe_copy)
        self.cycles_array = copy.deepcopy(self.cycles_array_copy)

    def bootstrap_atoms(self):
        if self.is_there_a_copy_of_atoms is False:
            self.save_copy_of_atoms()
            print(
                "[Warning] : I just saved a copy of the atom dataframe because you will destruct your original dataframe."
            )
        new_atoms = []
        for n in range(self.n_cycles):
            cycle = random.choice(self.cycles_array_copy)
            df = copy.deepcopy(
                self.atoms_dataframe_copy[self.atoms_dataframe_copy["Cycle"] == cycle]
            )
            df["Cycle"] = n * np.ones(len(df))
            new_atoms.append(df)
        self.atoms = pd.concat(new_atoms)
        self.cycles_array = self.atoms["Cycle"].unique()
        if len(self.cycles_array) != self.n_cycles:
            print("WWWWWWWAAAAAAAAAAAAA something went wrong it is weird.")

    def define_variable1(self, **kwargs):
        self.var1 = Variable(**kwargs)

    def define_variable2(self, **kwargs):
        self.var2 = Variable(**kwargs)

    def get_atoms_in_box(self, df, box):
        """
        Retourne un dataframe avec les positions de atomes à l'intérieur de la boîte.

        Parameters
        ----------
        df : dataframe d'atomes
        box : dictionnaire, du type {"Vx": {"size": 10, "position": 0}}. Il faut que les entrées du dictionnaire matchent le nom des colonnes du dataframe soit Vx, Vy, Vz et Cycle. Update depuis mai 2023, les entrées du dictionnaire de boite peuvent être 'range' (array) ou bien minimum et maximum.

        Returns
        ----------
        df : dataframe avec les même colonnes dont les atomes sont tous dans la box.
        """
        for key, value in box.items():
            # Rappel : key est par ex "Vx" ; value est {"size":10, "position":0}
            if "range" in value:
                minimum = np.min(value["range"])
                maximum = np.max(value["range"])
            elif "position" in value and "size" in value:
                minimum = value["position"] - np.abs(value["size"]) / 2
                maximum = value["position"] + np.abs(value["size"]) / 2
            elif "minimum" in value and "maximum" in value:
                minimum = value["minimum"]
                maximum = value["maximum"]
            df = df[((df[key] >= minimum) & (df[key] < maximum))]
        return df

    def apply_ROI(self):
        """
        Modifie le dataframe "atoms" en appliquant la ROI. Cela permet d'alléger les données à traiter.
        Si la ROI est vide, la méthode ne fait rien. Le format de la ROI doit être {"Vx": {"max":120, "min":-120}}
        """
        if self.ROI:
            for key, entry in self.ROI.items():
                self.atoms = self.atoms[
                    (
                        (self.atoms[key] <= entry["max"])
                        & (self.atoms[key] > entry["min"])
                    )
                ]

    def apply_ROD(self):
        """
        Modifie le dataframe "atoms" en appliquant la region of desinterest i.e. en sélectionnant les atomes autre que ceux dans la ROD. Si la ROF est vide, la méthode ne fait rien.
        """
        if self.ROD:
            for key, entry in self.ROD.items():
                self.atoms = self.atoms[
                    (
                        (self.atoms[key] > entry["max"])
                        | (self.atoms[key] < entry["min"])
                    )
                ]

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

    def counts_atoms_in_boxes_one_variable(self, df, var, box, column_name="N_1"):
        """
        Prend en argument un dataframe d'atomes, une variable et une boîte. Pour chaque "value" de  "variable", elle redéfinie la taille/position de la boîte et récupère le nombre d'atomes dans la boîte à chaque cycle. Elle renvoie un dataframe de 3 colonnes : une avec les cycles de nom "Colonne", une avec le nombre d'atome au cycle donné (de nom column_name) et une avec la valeur de la position/taille de la boîte (de nom var.name)

        Parameters
        ----------
        df : pandas dataframe, dataframe avec les atomes
        var : Variable
        box : dictionnaire, position et taille de la boîte sur 1 à 4 axes.
            Exemple {"Vx": {"size": 10, "position": 0}}

        Returns
        ----------
        result : pandas dataframe
            dataframe avec 3 colonnes :
                "Cycle" : avec le numéro du cycle
                column_name (ex : "N1") : avec le nombre d'atome dans la boîte
                var.name (ex : "ΔVx"): avec la valeur de la position/taille de la boîte.
        """
        # On parcourt les différentes valeurs de la Variable var (tqdm --> waiting bar)
        # for i in tqdm(range(var.n_step), desc="Gathering {}".format(var.name)):
        for i in range(var.n_step):
            # On change la boîte selon la i-ème valeur de var
            box[var.axe][var.type] = var.get_value_i(i)
            # On récupère le dataframe avec le nombre d'atome dans la boîte à chaque cycle
            dataframe = self.obtain_number_of_atoms_per_cycle_in_box(
                df, box, column_name=column_name
            )
            dataframe[var.name] = var.get_value_i(i)
            # On ajoute la colonne avec la valeur de variable
            if i == 0:
                result = dataframe
            else:
                result = pd.concat([result, dataframe])
        # les index se répètent : c'est bof : je les réinitialise.
        result.reset_index(drop=True)
        return result

    ###########################################################
    ############### FONCTION D'AFFICHAGE  #####################
    ###########################################################

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
        plt.figure()
        X_list = self.atoms[nameX].to_numpy()
        Y_list = self.atoms[nameY].to_numpy()
        hist_values, X_values, Y_values, _ = plt.hist2d(
            X_list, Y_list, bins=[x_bins, y_bins], cmap=plt.cm.Blues
        )
        if show_boxes:

            def draw_box(cX, σX, cY, σY, color="orange", label="box"):
                plt.plot(
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

            ### On affiche 3x les boites
            # on affiche la boite 1
            posX = self.boxes["1"][nameX]["position"]
            sizeX = 3 * self.boxes["1"][nameX]["size"]
            posY = self.boxes["1"][nameY]["position"]
            sizeY = 3 * self.boxes["1"][nameY]["size"]

            # et le boite 2
            posX = self.boxes["2"][nameX]["position"]
            sizeX = 3 * self.boxes["2"][nameX]["size"]
            posY = self.boxes["2"][nameY]["position"]
            sizeY = 3 * self.boxes["2"][nameY]["size"]
        plt.colorbar()
        plt.legend()
        plt.title(title)
        plt.xlabel(nameX)
        plt.ylabel(nameY)
        if save:
            plt.savefig(save)
        if show_plot is True:
            plt.show()
        else:
            plt.close()
        return (hist_values, X_values, Y_values)

    def show_2D_plot(self, x, y, z="g^2", vmin=0.8, vmax=2):
        """
        Show a 2D plot of some computed values.

        Parameters
        ----------
        x : string, colonne du dataframe result à plotter
        y : string, idem
        z : string, idem
        vmin : float, valeur minimale pour la colormap
        vmax : float, valeur maximale pour la colormap
        """
        df_pivoted_correlations = self.result.pivot(index=x, columns=y, values=z)
        fig, ax = plt.subplots(figsize=(13, 10))
        sns.heatmap(
            df_pivoted_correlations,
            cmap="YlGnBu",
            ax=ax,
            vmin=vmin,
            vmax=vmax,
        )
        ax.invert_yaxis()
        plt.title(z)
        plt.show()

    ###########################################################
    ######## GESTION DES CAS DE FIGURE AVANT LE CALCUL  #######
    ###########################################################

    def compute_correlations(self):
        """
        Cette fonction gère les différents cas de figure de scan.
        """
        # Cas 1 : on ne scanne aucun paramètre, on veut juste la corrélation entre deux boites.
        if (self.var1 == None) and (self.var2 == None):
            atoms_box1 = self.obtain_number_of_atoms_per_cycle_in_box(
                self.atoms, self.boxes["1"], column_name="N_1"
            )
            atoms_box2 = self.obtain_number_of_atoms_per_cycle_in_box(
                self.atoms, self.boxes["2"], column_name="N_2"
            )
            total_atoms = self.merge_dataframe_on_cycles(atoms_box1, atoms_box2)
            corr_names = self.quantity_of_interest()
            corr_values = self.quantity_of_interest(total_atoms)
            self.result = pd.DataFrame([corr_values], columns=corr_names)
        # Cas 2 : on ne scanne qu'un seul paramètre : dans ce cas, on définit la seconde variable scannée avec un seul paramètre en appelant define_new_variable_with_one_value puis on réappelle la méthode compute_correlations pour être dans le cas 3
        elif (self.var1 == None) and (self.var2 != None):
            # dans ce cas on va définir la variable 1
            # print("I defined myself variable1")
            self.define_new_variable_with_one_value(self.var2, 1)
            self.compute_correlations()
        elif (self.var1 != None) and (self.var2 == None):
            # print("I defined myself variable2")
            self.define_new_variable_with_one_value(self.var1, 2)
            self.compute_correlations()
        # Cas 3 : on scanne des paramètres appartenant à deux boîtes différentes
        elif self.var1.box != self.var2.box:
            self.compute_correlations_different_box_scanned()
        # Cas 4 : on scanne des paramètres appartenant à la même boîte
        else:
            self.result = 0
            print("This is not implemented yet.")

    def compute_correlations_different_box_scanned(self):
        """
        Méthode pour calcul des corrélations lorsque var1 et var2 (les paramètres scannés) correspondent à deux boites différentes.
        """
        # %#% STEP 1 : on récupère le nombre d'atome dans les boites associées à var1 et var2. Disons que intuivement, var1 corresponde à la boîte 1 et var2 à la boîte 2 mais ce n'est pas nécessaire dans le code.
        # --> Start with var1
        # On ne récupère que les atomes présents dans la boîte selon les deux axes qui ne varient pas. Par exemple, si on est en train de faire varier la position de la boite selon Vx, on récupère les atomes qui vérifient déjà les bonnes conditions selon Vy et Vz pour alléger les calculs.
        box = self.boxes[self.var1.box].copy()
        # On enlève l'axe concerné par le scan
        posi_and_size = box.pop(self.var1.axe)
        # l'axe concerné par le scan est donc
        scanned_box = {self.var1.axe: posi_and_size}
        # On a donc deux boîtes maintenant : une avec deux axes (box) qui ne sont pas modifié à chaque position/taille de boîte et une autre (scanned_box) avec un seul axe qui correspond à l'axe scanné var1.axe
        df_atoms_var1 = self.get_atoms_in_box(self.atoms, box)
        # Result var1 est un dataframe avec le nombre d'atomes dans la boîte à chaque cycle pour les différentes positions de la boîte. Il a donc 3 colonnes dont les entêtes sont "Cycle", "N_i" avec i = 1 ou 2 selon la boîte et le nom de la variables scannée par exemple "ΔVx". Voir la documentation de counts_atoms_in_boxes_one_variable pour plus de détails.
        result_var1 = self.counts_atoms_in_boxes_one_variable(
            df_atoms_var1, self.var1, scanned_box, column_name="N_" + self.var1.box
        )
        ## check 21 of may:
        column_name = "N_" + self.var1.box
        df = result_var1[result_var1[column_name] < 0]
        if len(df) > 1:
            print(df)
        # --> Do the same with var2
        box = self.boxes[self.var2.box].copy()
        posi_and_size = box.pop(self.var2.axe)
        scanned_box = {self.var2.axe: posi_and_size}
        df_atoms_var2 = self.get_atoms_in_box(self.atoms, box)
        result_var2 = self.counts_atoms_in_boxes_one_variable(
            df_atoms_var2, self.var2, scanned_box, column_name="N_" + self.var2.box
        )
        column_name = "N_" + self.var2.box
        df = result_var2[result_var2[column_name] < 0]
        if len(df) > 1:
            print(df)
        # %#% STEP2
        # On construit le dataframe total, qui initialement contient 5 colonnes : Cycle le cycle, N_1 et N_2 nombre d'atomes dans la boîte 1 et 2, self.var1 et self.var2 la position/taille des boîtes lors du scan. Le nombre de lignes de total est dont Nombre_de_cycles x Nombre_de_différentes_var1 x Nombre_de_différentes_var2.
        total = pd.merge(result_var1, result_var2, on="Cycle")

        # %#% STEP3 : computes quantites of interest
        self.compute_result(total)

    def compute_opposite_momenta_correlations(self):
        """
        Calcule des corrélation dans self.result lorsqu'on ne veut scanner qu'une variable et faire des corrélations entre k et -k. Les paramètres qui sont bougés doivent être définis dans var1.
        """
        box = self.var1.box
        if box == "1":
            box = "2"
        else:
            box = "1"
        axe = self.var1.axe
        type = self.var1.type
        name = self.var1.name
        values = self.var1.values
        self.define_variable2(
            box=box, axe=axe, type=type, name="-" + name, values=-1 * values
        )

        # On fait maintenant la même chose que dans compute_correlations_different_box_scanned() jusqu'à step 2
        #### STEP 1 : on récupère le nombre d'atome dans les boites associées à var1 et var2. Disons que intuivement, var1 corresponde à la boîte 1 et var2 à la boîte 2 mais ce n'est pas nécessaire dans le code.
        # --> Start with var1
        # On ne récupère que les atomes présents dans la boîte selon les deux axes qui ne varient pas. Par exemple, si on est en train de faire varier la position de la boite selon Vx, on récupère les atomes qui vérifient déjà les bonnes conditions selon Vy et Vz pour alléger les calculs.
        box = self.boxes[self.var1.box].copy()
        # On enlève l'axe concerné par le scan
        posi_and_size = box.pop(self.var1.axe)
        # l'axe concerné par le scan est donc
        scanned_box = {self.var1.axe: posi_and_size}
        # On a donc deux boîtes maintenant : une avec deux axes (box) qui ne sont pas modifié à chaque position/taille de boîte et une autre (scanned_box) avec un seul axe qui correspond à l'axe scanné var1.axe
        df_atoms_var1 = self.get_atoms_in_box(self.atoms, box)
        # Result var1 est un dataframe avec le nombre d'atomes dans la boîte à chaque cycle pour les différentes positions de la boîte. Il a donc 3 colonnes dont les entêtes sont "Cycle", "N_i" avec i = 1 ou 2 selon la boîte et le nom de la variables scannée par exemple "ΔVx". Voir la documentation de counts_atoms_in_boxes_one_variable pour plus de détails.
        result_var1 = self.counts_atoms_in_boxes_one_variable(
            df_atoms_var1, self.var1, scanned_box, column_name="N_" + self.var1.box
        )
        # --> Do the same with var2
        box = self.boxes[self.var2.box].copy()
        posi_and_size = box.pop(self.var2.axe)
        scanned_box = {self.var2.axe: posi_and_size}
        df_atoms_var2 = self.get_atoms_in_box(self.atoms, box)
        result_var2 = self.counts_atoms_in_boxes_one_variable(
            df_atoms_var2, self.var2, scanned_box, column_name="N_" + self.var2.box
        )
        #### STEP2
        result_var1["abs(varname)"] = np.abs(result_var1[self.var1.name])
        result_var2["abs(varname)"] = np.abs(result_var2[self.var2.name])
        # On fusionne maintenant les dataframe sur deux clés et non plus sur une seule.
        total = pd.merge(result_var1, result_var2)

        #### STEP3 : computes quantites of interest
        self.compute_result(total)

    def define_new_variable_with_one_value(self, var, new_var_number):
        """Cette fonction génère la variable 1 ou 2 (correspondant à var_number) en copiant les paramètres de var. Par exemple, on scanne la variable 1, la taille de boîtes selon Vz, cette fonction va définir la variable 2 comme scannant la taille des boites dans la boite opposée à la variable 1 avec un seul paramètre, celui définit dans box. Cela permet d'utiliser systématiquement la fonction compute_correlations_different_box_scanned même lorsuq'on a une seule boîte définie.

        Parameters
        ----------
        var : objet de la classe Variable

        var_number : entier 1 ou 2
            Si c'est 1, on va définir une nouvelle variable 1, si c'est 2, on va définir la nouvelle variable 2.
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
        Prend en argument un dataframe dont les colonnes sont le nombre d'atomes dans les boîtes 1 et 2 à chaque cycle1 Retourne une liste de flottants avec les différentes valeurs d'intérêt. Si dataframe = "", renvoie une liste de strings avec le nom de chaque quantité calculée.
        Pour rajouter des valeurs d'intérêt à calculer, il suffit de rajouter dans la liste column_names le nom de la variable et de l'ajouter dans result.

        Parameters
        ----------
        dataframe : pandas dataframe avec (au moins) 3 colonnes
                "Cycle" le numéro de cycle
                "N_1", avec le nombre d'atomes dans la boîte 1
                "N_2", avec le nombre d'atomes dans la boîte 2

        Returns
        ----------
        column_names (si aucun dataframe n'est donné) : liste de string avec le nom des valeurs calculées
        computed_quantities (si un dataframe est donné) : liste de flottants avec la quantité d'intérêt calculée.
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
        Calcul le dataframe de corrélation result utilisant total. Total est un dataframe avec pour chaque cycle et chaque position/taille de boites le nombre d'atomes dans la boite 1 et dans la boite 2.

        /!\ /!\ Ne fonctionne que si on a deux variables --> à faire pour un scan 1D ? ou bien on garde la façon naïve car ce n'est pas long.

        Parameters
        ----------
        total : pandas dataframe de 5 colonnes : 'Cycle' le cycle, 'N_1' et 'N_2' nombre d'atomes dans la boîte 1 et 2, self.var1 et self.var2 la position/taille des boîtes lors du scan. Le nombre de lignes de total est dont Nombre_de_cycles x Nombre_de_différentes_var1 x Nombre_de_différentes_var2.


        Built
        ----------
        self.result : pandas dataframe de Nombre_de_différentes_var1 x Nombre_de_différentes_var2 lignes. Ses différentes colonnes sont
            var1.name et var2.name : le nom des variables scannées
            "N_1", "N_2" :  moyenne sur les cycles du nombre d'atomes dans la boîte 1, 2
            "N_1+N_2" : somme des populations moyenne
            "N_1-N_2" : différence des populations moyenne
            "(N_1-N_2)^2" : moyenne (sur les cycles) des écart quadratiques
            "N_1*N_2" : moyenne du produit du nombre d'atomes
            "variance" : variance, il s'agit de la quantité <(N1 - N2)^2> / <N1-N2>^2 où la moyenne est prise sur les cycles.
            "normalized variance" : variance normalisée, variance divisée par la somme des populations  : variance / <N1 + N2>
            "<N_1>*<N_2>" : produit de la moyenne du nombre d'atomes
        """
        total["N_1*N_2"] = total["N_1"] * total["N_2"]
        column_name = "N_1*N_2"
        total["N_1**2"] = total["N_1"] ** 2
        total["N_2**2"] = total["N_2"] ** 2
        total["N_1-N_2"] = total["N_1"] - total["N_2"]
        total["(N_1-N_2)^2"] = (total["N_1"] - total["N_2"]) ** 2
        total["N_1+N_2"] = total["N_1"] + total["N_2"]
        total["(N_1-N_2)^4"] = (total["N_1"] - total["N_2"]) ** 4
        self.total = total
        # on moyenne les données sur les cycles (on les groupe donc par différentes valeurs de var1 et var2)
        self.result = total.groupby(
            [self.var1.name, self.var2.name], as_index=False
        ).mean()
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

        # ---------------
        # Variance
        # ---------------
        self.result["variance"] = (
            self.result["(N_1-N_2)^2"] - (self.result["N_1-N_2"]) ** 2
        )

        self.result["normalized variance"] = self.result["variance"] / (
            self.result["N_1"] + self.result["N_2"]
        )

        # ---------------
        # Calculs de g^2
        # ---------------
        self.result["g^2"] = self.result["N_1*N_2"] / (
            self.result["N_1"] * self.result["N_2"]
        )

        # Cauchy Schwarz
        self.result["C-S"] = self.result["N_1*N_2"] / (
            np.sqrt(
                (self.result["N_1**2"] - self.result["N_1"])
                * (self.result["N_2**2"] - self.result["N_2"])
            )
        )
        self.result["g^2(k1,k1)"] = (
            self.result["N_1**2"] - self.result["N_1"]
        ) / self.result["N_1"] ** 2
        self.result["g^2(k2,k2)"] = (
            self.result["N_2**2"] - self.result["N_2"]
        ) / self.result["N_2"] ** 2

        # on enlève le shot noise si cela est demandé par l'utilisateur.
        if self.remove_shot_noise:
            if self.var1.type == "position":
                VJ = self.var1.axe  # I assume that we scanned the same axe
                local_condition = (
                    self.result[self.var1.name] == self.result[self.var2.name]
                )
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
                    if len(local_condition) > 0:
                        print(
                            "[WARNING] Shot Noise has not been taken off weird because boxes are do not have the same center. "
                        )
                elif (
                    (self.boxes["1"]["Vz"]["size"] != self.boxes["2"]["Vz"]["size"])
                    or (self.boxes["1"]["Vy"]["size"] != self.boxes["2"]["Vy"]["size"])
                    or (self.boxes["1"]["Vx"]["size"] != self.boxes["2"]["Vx"]["size"])
                ):
                    if len(local_condition) > 0:
                        print(
                            "[WARNING] Shot Noise has not been taken off weird because boxes do not have the same size. Please be carefull when delaing with local correlations !"
                        )
                else:
                    self.result.loc[local_condition, "g^2"] = (
                        self.result["N_1*N_2"] - self.result["N_1"]
                    ) / (self.result["N_1"] * self.result["N_2"])
                    self.result.loc[local_condition, "N_1*N_2"] = (
                        self.result["N_1*N_2"] - self.result["N_1"]
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

        # ---------------
        # Déviations standards
        # ---------------

        self.result["N_1 std"] = np.sqrt(self.result["N_1"])
        self.result["N_2 std"] = np.sqrt(self.result["N_2"])
        self.result["N_1 rel"] = self.result["N_1 std"] / self.result["N_1"]
        self.result["N_2 rel"] = self.result["N_2 std"] / self.result["N_2"]
        self.result["N_1*N_2 std"] = self.result["N_1*N_2"] * np.sqrt(
            self.result["N_1 rel"] ** 2 + self.result["N_2 rel"] ** 2
        )
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

        self.result["N_1 error"] = self.result["N_1 std"] / np.sqrt(self.n_cycles)
        self.result["N_2 error"] = self.result["N_2 std"] / np.sqrt(self.n_cycles)
        self.result["N_1*N_2 error"] = self.result["N_1*N_2 std"] / np.sqrt(
            self.n_cycles
        )
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
        # On exprime les densités avec des unités
        # ---------------
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
        # ---------------
        # On rajoute la différence et la moyenne des var1 et var2
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

        # print("Computation is done.")

    def get_atoms_distribution(
        self, nbMax, nbPt, posZ, sizeZ, posX, sizeX, posY, sizeY
    ):
        """
        Permet de tracer la distribution du nombre d'atomes moyenne sur tous les cycles, soit dans une boîte soit dans une moyenne de boîtes
        Parameters
        ----------
        posZ, sizeZ, posX, sizeX, posY, sizeY : paramètres de la boîte initiale
        nbMax : nombre d'atomes simultanés détectés à considérer
        nbPt : nombre de boîtes sur lesquelles on moyenne (la 1e boîte est celle centrée sur PosZ et de taille sizeZ, puis les autres boîtes sont les voisines selon Vz dans le sens des Vz croissants et de taille sizeZ

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
        self, nbMax, nbPt, posZ1, posZ2, sizeZ, posX, sizeX, posY, sizeY, show=True
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

        # fonction auxiliaire pour construire les boîtes
        def return_boxes(i, posZ1, posZ2, posX, posY, sizeZ, sizeX, sizeY):
            assert nbPt % 2 == 1
            if i == 0:
                box1_proba = {
                    "Vz": {"position": posZ1, "size": sizeZ},
                    "Vy": {"position": posY, "size": sizeY},
                    "Vx": {"position": posX, "size": sizeX},
                }
                box2_proba = {
                    "Vz": {"position": posZ2, "size": sizeZ},
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

    def return_dictionary_correlation_property(self) -> dict:
        """Return a dictionay with all the parameter of the simulation.

        Returns
        -------
        dict
            dictionary with all parameters of the correlation.
        """
        from flatten_dict import flatten, reducers

        dictionary = {}
        for key, value in self.__dict__.items():
            if type(value) in [int, float, dict, bool]:
                dictionary[key] = value
            elif type(value) == Variable:
                dictionary[key] = {}
                val_dic = copy.deepcopy(value.__dict__)
                del val_dic["values"]
                dictionary[key] = val_dic
        dictionary = flatten(dictionary, reducer=reducers.make_reducer(delimiter=" | "))
        return dictionary

    def return_pandas_dataframe_correlation_properties(self) -> pd.DataFrame:
        dictionary = self.return_dictionary_correlation_property()
        df = pd.DataFrame(data=[dictionary.values()], columns=dictionary.keys())
        return df


class Variable:
    """
    Définie la classe Variable.
    """

    def __init__(self, **kwargs):
        self.box = "1"  # le numéro de la boîte du paramètre scanné (1 ou 2)
        self.axe = "Vx"  # l'axe du paramètre scanné (Vx, Vy ou Vz)
        self.type = "size"  # le type : size ou position
        self.name = "ΔVx"  # son nom pour la colonne dans le dataframe
        self.min = 0  # sa valeur minimale, maximale et le nombre d'éléments
        self.max = 2
        self.step = 4
        self.values = []  # ou bien les valeurs qu'il prend (liste ou numpy array)
        self.round_decimal = 7
        self.__dict__.update(kwargs)
        if len(self.values) == 0:
            self.built_values()
        self.get_values_caracteristics()

    def built_values(self):
        self.values = np.arange(self.min, self.max, self.step)
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
        self.values = np.array(self.values)
        self.min = np.min(self.values)
        self.max = np.max(self.values)
        self.n_step = len(self.values)

    def get_value_i(self, i):
        """
        Renvoie la i-ième value de self.values
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
