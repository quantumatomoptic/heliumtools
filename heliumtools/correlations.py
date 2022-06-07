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
import copy


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

    def __init__(self, atoms, n_cycles, **kwargs):
        """
        Object initialization, sets parameters as the user defined, build the atoms dataframe and apply ROD and ROI.
        """
        self.atoms = atoms.copy()
        self.n_cycles = n_cycles
        self.bec_arrival_time = 308  # temps d'arrivée du BEC, en ms
        self.raman_kick = 42  # mm/s, kick Raman
        self.gravity = 9.81
        self.var1 = None
        self.var2 = None
        self.remove_shot_noise = True
        self.ROI = {}  # Region Of Interest
        self.ROD = {}  # Region Of Desinterest.
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
        self.__dict__.update(kwargs)
        self.boxes = self.boxes.copy()
        self.initial_boxes = copy.deepcopy(self.boxes)  # perform a deep copy
        # atoms : (mm,mm, ms) --> (mm/s, mm/s, mm/s)
        self.build_the_atoms_dataframe()
        self.apply_ROI()  # Keep only atoms in ROI
        self.apply_ROD()  # Take off atoms in Region Of Desinterest.
        print("Data are loaded")
        # Initialisation de self.result datframe avec toutes les corrélations
        # self.result = pd.DataFrame(
        #     np.zeros(1, len(self.quantity_of_interest())),
        #     columns=self.quantity_of_interest,
        # )

    def build_the_atoms_dataframe(self):
        """
        Cette méthode construit le dataframe contenant l'ensemble des positions des atomes : elle construit le dataframe self.atoms, dont les vitesses sont exprimées en mm/s à partir du dataframe initial.
        """
        self.atoms["X"] = 1000 * self.atoms["X"] / self.atoms["T"] + self.raman_kick
        self.atoms = self.atoms.rename(columns={"X": "Vx"})
        self.atoms["Y"] = 1000 * self.atoms["Y"] / self.atoms["T"]
        self.atoms = self.atoms.rename(columns={"Y": "Vy"})
        self.atoms["T"] = self.gravity * (self.atoms["T"] - self.bec_arrival_time)
        self.atoms = self.atoms.rename(columns={"T": "Vz"})

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
        box : dictionnaire, du type {"Vx": {"size": 10, "position": 0}}. Il faut que les entrées du dictionnaire matchent le nom des colonnes du dataframe soit Vx, Vy, Vz et Cycle.

        Returns
        ----------
        df : dataframe avec les même colonnes dont les atomes sont tous dans la box.
        """
        for key, value in box.items():
            # Rappel : key est par ex "Vx" ; value est {"size":10, "position":0}
            minimum = value["position"] - value["size"] / 2
            maximum = value["position"] + value["size"] / 2
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
        my_serie = df.value_counts(subset="Cycle")

        atoms_in_box = my_serie.to_frame()
        atoms_in_box = atoms_in_box.rename(columns={0: column_name})
        atoms_in_box.reset_index(inplace=True)
        # atoms_in_box is now a dataframe with two columns "Cycle" and "N_1" (or column_name). However, if there were no atom at cycle 34 in the box, this cycle does not appear inside atoms_in_box. In order to have the number of atoms in the box at each cycle, we must add 0 to those cycles which does not appear.
        # cycle_dataframe is just a dataframe with n_cycles : we use it to merge and add zeros to atoms_in_box
        cycle_dataframe = pd.DataFrame(
            np.arange(1, self.n_cycles + 1, 1), columns=["Cycle"]
        )
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
        for i in tqdm(range(var.n_step), desc="Gathering {}".format(var.name)):
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

    def compute_correlations_different_box_scanned(self):
        """
        /!\ Deprecated because it is too slow.... /!\

        Méthode pour calcul des corrélations lorsque var1 et var2 (les paramètres scannés) correspondent à deux boites différentes.
        """
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
        # On construit un gros dataframe de plusieurs colonnes, deux donnant les valeurs des positions/tailles des boîtes scannées et le reste étant des variables calculées par la méthode quantity_of_interest. Le nom des colonnes des variables calculées par quantity_of_interest est obtenu en appelant cette méthode sans argument.
        column_name = [self.var1.name, self.var2.name] + self.quantity_of_interest()
        column_values = np.zeros(
            (self.var1.n_step * self.var2.n_step, len(column_name))
        )
        self.result = pd.DataFrame(column_values, columns=column_name)
        step = 0
        with alive_bar(
            self.var1.n_step * self.var2.n_step, title="Computing correlations"
        ) as bar:
            for i1 in range(self.var1.n_step):
                for i2 in range(self.var2.n_step):
                    var1_value = self.var1.get_value_i(i1)
                    var2_value = self.var2.get_value_i(i2)
                    # On sélectionne les parties des dataframe où les variables ont la bonne valeur (cemme de l'étape i1, i2) avant de les merger et de les corréler.
                    var1dataframe = result_var1[
                        result_var1[self.var1.name] == var1_value
                    ]
                    var2dataframe = result_var2[
                        result_var2[self.var2.name] == var2_value
                    ]
                    dataframe = self.merge_dataframe_on_cycles(
                        var1dataframe, var2dataframe
                    )
                    line_result = [var1_value, var2_value] + self.quantity_of_interest(
                        dataframe
                    )
                    self.result.loc[self.result.index[step]] = line_result
                    step += 1
                    bar()

    def compute_correlations_one_box_scanned(self, var):
        """
        Calcule des corrélations self.result lorsqu'une seule boîte est scannée. La démarche est la même que celle de compute_correlations_different_box_scanned
        """
        # On ne récupère que les atomes présents dans la boîte selon les deux axes qui ne varient pas. Par exemple, si on est en train de faire varier la position de la boite selon Vx, on récupère les atomes qui vérifient déjà les bonnes conditions selon Vy et Vz pour alléger les calculs.
        box = self.boxes[var.box].copy()
        # On enlève l'axe concerné par le scan
        posi_and_size = box.pop(var.axe)
        # l'axe concerné par le scan est donc
        scanned_box = {var.axe: posi_and_size}
        # On a donc deux boîtes maintenant : une avec deux axes (box) qui ne sont pas modifié à chaque position/taille de boîte et une autre (scanned_box) avec un seul axe qui correspond à l'axe scanné var1.axe
        df_atoms_var = self.get_atoms_in_box(self.atoms, box)
        result_var = self.counts_atoms_in_boxes_one_variable(
            df_atoms_var, var, scanned_box, column_name="N_" + var.box
        )
        if var.box == "1":
            result_not_var = self.obtain_number_of_atoms_per_cycle_in_box(
                self.atoms, self.boxes["2"], column_name="N_2"
            )
        elif var.box == "2":
            result_not_var = self.obtain_number_of_atoms_per_cycle_in_box(
                self.atoms, self.boxes["1"], column_name="N_1"
            )

        column_name = [var.name] + self.quantity_of_interest()
        column_values = np.zeros((var.n_step * 1, len(column_name)))
        self.result = pd.DataFrame(column_values, columns=column_name)
        for i1 in range(var.n_step):
            var_value = var.get_value_i(i1)

            dataframe = self.merge_dataframe_on_cycles(
                result_var[result_var[var.name] == var_value], result_not_var
            )
            line_result = [var_value] + self.quantity_of_interest(dataframe)
            self.result.loc[self.result.index[i1]] = line_result

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
        # Cas 2 : on ne scanne qu'un seul paramètre
        elif (self.var1 == None) and (self.var2 != None):
            self.compute_correlations_one_box_scanned(self.var2)
        elif (self.var1 != None) and (self.var2 == None):
            self.compute_correlations_one_box_scanned(self.var1)
        # Cas 3 : on scanne des paramètres appartenant à deux boîtes différentes
        elif self.var1.box != self.var2.box:
            self.compute_correlations_different_box_scanned_fast()
        # Cas 4 : on scanne des paramètres appartenant à la même boîte
        else:
            self.result = 0
            print("This is not implemented yet.")

        self.boxes = self.initial_boxes

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
        total["N_1-N_2"] = total["N_1"] - total["N_2"]
        total["(N_1-N_2)^2"] = (total["N_1"] - total["N_2"]) ** 2
        total["N_1+N_2"] = total["N_1"] + total["N_2"]
        total["(N_1-N_2)^4"] = (total["N_1"] - total["N_2"]) ** 4
        self.total = total
        # on moyenne les données sur les cycles (on les groupe donc par différentes valeurs de var1 et var2)
        self.result = total.groupby(
            [self.var1.name, self.var2.name], as_index=False
        ).mean()
        # On fait ensuite une petite manipulation pour calculer l'erreur sur la variance.
        # L'idée est de définir la variable (N_1-N_2)^2-moy(N_1-N_2)^2 puis de dire à l'ordinateur de calculer sa variance tout seul pour éviter de mettre la formule très longue et compliquée (wiki du 2 juin 2022). On reconstruit un dataframe avec les numéros de cycles
        df1 = self.result[[self.var1.name, self.var2.name, "N_1-N_2"]]
        df2 = pd.DataFrame({"Cycle": np.linspace(1, self.n_cycles, self.n_cycles)})
        new_df = pd.merge(df1, df2, how="cross")
        new_df = new_df[["Cycle", self.var1.name, self.var2.name, "N_1-N_2"]]
        new_df.columns = new_df.columns.str.replace("N_1-N_2", "mean(N_1-N_2)")
        # new_df est donc un dataframe avec 4 colonnes : une kz1, une kz2, une avec le cycle et une avec la moyenne (N1-N2). Bien entendu, à chaque cycle la moyenne N1-N2 est la même. NB : kz1 est de façon générale var1.name mais souvent kz1.
        # On veut ajouter au dataframe total la colonne moy(N1-N2).
        total = pd.merge(total, new_df)
        total["(N_1-N_2)^2-mean(N_1-N_2)^2"] = (
            total["(N_1-N_2)^2"] - total["mean(N_1-N_2)"] ** 2
        )

        error = total.groupby([self.var1.name, self.var2.name], as_index=False).std()

        print("Total dataframe is summed already")

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

        # on enlève le shot noise si cela est demandé par l'utilisateur.
        if self.remove_shot_noise:
            local_condition = self.result[self.var1.name] == self.result[self.var2.name]
            self.result.loc[local_condition, "g^2"] = (
                self.result["N_1*N_2"] - self.result["N_1"]
            ) / (self.result["N_1"] * self.result["N_2"])

        # ---------------
        # Déviations standards
        # ---------------
        self.result["N_1 std"] = error["N_1"]
        self.result["N_2 std"] = error["N_2"]
        self.result["N_1*N_2 std"] = error["N_1*N_2"]
        self.result["N_1-N_2 std"] = error["N_1-N_2"]
        self.result["(N_1-N_2)^2 std"] = error["(N_1-N_2)^2"]
        self.result["N_1+N_2 std"] = error["N_1+N_2"]
        self.result["variance std"] = error["(N_1-N_2)^2-mean(N_1-N_2)^2"]

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

        print("Computation is done.")

    def compute_correlations_different_box_scanned_fast(self):
        """
        Méthode pour calcul des corrélations lorsque var1 et var2 (les paramètres scannés) correspondent à deux boites différentes.
        """
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
        # On construit le dataframe total, qui initialement contient 5 colonnes : Cycle le cycle, N_1 et N_2 nombre d'atomes dans la boîte 1 et 2, self.var1 et self.var2 la position/taille des boîtes lors du scan. Le nombre de lignes de total est dont Nombre_de_cycles x Nombre_de_différentes_var1 x Nombre_de_différentes_var2.
        total = pd.merge(result_var1, result_var2, on="Cycle")

        #### STEP3 : computes quantites of interest
        self.compute_result(total)


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
        self.values = None  # ou bien les valeurs qu'il prend (liste ou numpy array)
        self.round_decimal = 7
        self.__dict__.update(kwargs)
        if self.values == None:
            self.built_values()
        self.get_values_caracteristics()

    def built_values(self):
        self.values = np.arange(self.min, self.max, self.step)
        self.values = np.round(self.values, self.round_decimal)

    def get_values_caracteristics(self):
        self.values = np.array(self.values)
        self.min = np.min(self.values)
        self.max = np.min(self.values)
        self.n_step = len(self.values)

    def get_value_i(self, i):
        """
        Renvoie la i-ième value de self.values
        """
        return self.values[i]


import os, glob, math, random
from pathlib import Path


def load_data_for_correlation(folder, maximum_number_of_files=1e8):
    """
    Charge les .atoms d'une séquence donnée
    """
    ext = [".atoms"]
    selected_files = sorted(
        [
            path.as_posix()
            for path in filter(lambda path: path.suffix in ext, folder.glob("*"))
        ]
    )

    N_files = min(len(selected_files), maximum_number_of_files)
    Xc, Yc, Tc, Cyclec = [], [], [], []
    for i in range(N_files):
        path = selected_files[i]
        X, Y, T = get_data(path)
        cycle = np.ones(len(X)) * (1 + i)
        Xc, Yc, Tc, Cyclec = combine(Xc, Yc, Tc, Cyclec, X, Y, T, cycle)
    # Je crée maintenant un dataframe
    data = np.transpose(np.array([Cyclec, Xc, Yc, Tc]))
    df = pd.DataFrame(data, columns=["Cycle", "X", "Y", "T"])
    return N_files, df


def combine(Xc, Yc, Tc, Cyclec, X, Y, T, cycle):
    """
        Combine les positions et temps d'arrivée des atomes plus le cycle
    correspondant aux fichier.
    """
    Xc = np.concatenate([Xc, X])
    Yc = np.concatenate([Yc, Y])
    Tc = np.concatenate([Tc, T])
    Cyclec = np.concatenate([Cyclec, cycle])
    return (Xc, Yc, Tc, Cyclec)


def get_data(path):
    v_perp_x = 1.02  # mm/ns
    v_perp_y = 1.13  # mm/ns
    time_resolution = 1.2e-10
    # time_to_pos = 2 * 0.98e-9

    atoms_file = np.fromfile(path, dtype="uint64")

    atoms = atoms_file * time_resolution

    events_list = atoms.reshape(int(len(atoms) / 4), 4).T

    Xmcp = 0.5 * v_perp_x * 1e9 * (events_list[1] - events_list[0])
    Ymcp = 0.5 * v_perp_y * 1e9 * (events_list[3] - events_list[2])

    X = (Xmcp + Ymcp) / np.sqrt(2)
    Y = (Ymcp - Xmcp) / np.sqrt(2)
    T = (events_list[0] + events_list[1] + events_list[2] + events_list[3]) / 4

    T = T * 1e3

    return (X, Y, T)
