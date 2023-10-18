#!/usr/bin/env python
# -*- mode:Python; coding: utf-8 -*-

# ----------------------------------
# Created on the 5-4-2023 by Victor
#
# Developped by Victor, ...
#
# Last (big) change on the 14 of April by Victor : replacing numpy istogramdd by torch, wich is 10 times faster than [5]
#
# Copyright (c) 2023 - Helium1@LCF
# ----------------------------------
#
"""
Content of correlations2.py
-----------------------------

Please document your code ;-).

Bibliography for the futur
A fast histogram was implemented better than numpy in [1]. However, the problem is that it does not take into account 3D datas.

---- Améliorations

Je suis tombé sur [4] où une personne (du CERN) n'est pas non plus satisfaite des histogram à N dimensions de numpy. Je le teste sur ma machine et le programme marche 10 fois plus vite. Cependant, cela ne suffit pas pour que je puisse faire tourner des corrélations sur mon ordinateur.

[1] https://pypi.org/project/fast-histogram/
[2] https://github.com/vaexio/vaex
[3] https://towardsdatascience.com/beyond-pandas-spark-dask-vaex-and-other-big-data-technologies-battling-head-to-head-a453a1f8cc13
[4] https://github.com/pytorch/pytorch/issues/29209
[5] https://pytorch.org/docs/1.11/generated/torch.histogramdd.html
"""
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# import snoop
from random import choice, shuffle, randint
from tqdm import tqdm, trange
import copy
import time
import torch
from .misc.gather_data import apply_ROI, apply_ROD
from .data_builder import DataBuilder


class CorrelationHe2Style(DataBuilder):
    """
    CorrelationHe2Style class that inherits from the DatBuilder class
    """

    def __init__(self, atoms, **kwargs):
        """
        Object initialization, sets parameters as the user defined, build the atoms dataframe and apply ROD and ROI.
        """
        super().__init__(atoms, **kwargs)
        self.only_one_beam = False
        self.computer_performance = 1
        self.precision = 100
        self.axis = ["Vx", "Vy", "Vz"]
        self.voxel_size = {"Vx": 0.1, "Vy": 0.1, "Vz": 0.1}
        self.voxel_numbers = {"Vx": 11, "Vy": 11, "Vz": 11}
        self.beams = {
            "A": {
                "Vx": {"size": 20, "position": 0},
                "Vy": {"size": 20, "position": 0},
                "Vz": {"size": 5, "position": 25},
            },
            "B": {
                "Vx": {"size": 20, "position": 0},
                "Vy": {"size": 20, "position": 0},
                "Vz": {"size": 5, "position": -25},
            },
        }
        self.denominator_ratio = 0.5
        self.__dict__.update(kwargs)

    ############################################
    ## BASIC FUNCTIONS CALLED BY THE PROGRAMM ##
    ############################################

    def update_atoms_in_beams(self):
        """
        This method updates atoms in beam A and B. It is this function that multiplies transform velocities from float to int so that calculations are faster.
        """
        self.atoms.reset_index(drop=True, inplace=True)
        self.atoms["index"] = np.arange(0, len(self.atoms))
        self.atomsA = copy.deepcopy(apply_ROI(self.atoms, self.beams["A"]))
        self.atomsB = copy.deepcopy(apply_ROI(self.atoms, self.beams["B"]))
        cols_to_remove = [
            col
            for col in self.atomsA.columns
            if col not in self.axis + ["Cycle", "index"]
        ]
        self.atomsA.drop(cols_to_remove, axis=1, inplace=True)
        cols_to_remove = [
            col
            for col in self.atomsB.columns
            if col not in self.axis + ["Cycle", "index"]
        ]
        self.atomsB.drop(cols_to_remove, axis=1, inplace=True)
        # we change the type of atoms data from float64 to integer (int16, int32) as small as possible.
        for Vj in self.axis:
            self.atomsA[Vj] = self.atomsA[Vj].astype(np.float32)
            self.atomsB[Vj] = self.atomsB[Vj].astype(np.float32)
        for col in ["Cycle", "index"]:
            self.atomsA[col] = self.atomsA[col].astype(np.int32)
            self.atomsB[col] = self.atomsB[col].astype(np.int32)

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

    def initialize_voxel_map_properties(self):
        """_summary_"""
        self.voxel_map_size = tuple(
            [int(self.voxel_numbers[axis]) for axis in self.axis]
        )
        # self.voxel_map = np.zeros(self.voxel_map_size)
        # self.voxel_size = {"Vx": 0.1, "Vy": 0.1, "Vz": 0.1}
        self.voxel_map_range = tuple(
            [
                (
                    -self.voxel_numbers[axis] * self.voxel_size[axis] / 2,
                    self.voxel_numbers[axis] * self.voxel_size[axis] / 2,
                )
                for axis in self.axis
            ]
        )
        voxel_centers = []
        for axis in self.axis:
            mini = (
                -self.voxel_numbers[axis] * self.voxel_size[axis] / 2
                + self.voxel_size[axis] / 2
            )
            maxi = +self.voxel_numbers[axis] * self.voxel_size[axis] / 2
            voxel_centers.append(np.arange(mini, maxi, step=self.voxel_size[axis]))
            if len(voxel_centers[-1]) != self.voxel_numbers[axis]:
                print(len(voxel_centers[-1]))
                raise Exception(
                    "Something strange appended in the voxel map initialization."
                )
        # create the result dataframe

        data = np.array(
            [
                [x, y, z]
                for x in voxel_centers[0]
                for y in voxel_centers[1]
                for z in voxel_centers[2]
            ]
        )
        self.result = pd.DataFrame(data=data, columns=self.axis)
        if self.only_one_beam is True:
            for column in [
                "G2AA",
                "G2AA denominator",
            ]:
                self.result[column] = np.zeros(len(self.result))
        else:
            for column in [
                "G2AA",
                "G2BB",
                "G2AB",
                "G2AA denominator",
                "G2BB denominator",
                "G2AB denominator",
                # "G2AB denominator1",
                # "G2AB denominator2",
            ]:
                self.result[column] = np.zeros(len(self.result))

    def get_G2(
        self, atX: pd.DataFrame, atY: pd.DataFrame, local=True, numerator=True
    ) -> pd.DataFrame:
        """Function that compute the 3D velocity difference between all atoms in the crossed atX x atY dataframe (atoms in beam X, Y being A or B). If local is True, it computes the difference while if local is False, it computes the sum.


        Detailed explanation with an example :
        atX = DataFrame("Vx":[1,2] ; "Vy":[9,8] ; "Vz":[0,0], "index": [1018, 1019], "Cycle":[112,112] )
        atY = atX
        We want to compute the local correlations. First step is to create the cross product of this DataFrame using merge. We create the following dataframe (_x and_y suffix are added by pandas) :
        atXY = "Vx_x", "Vy_x", "Vz_x", "index_x", "Cycle", "Vx_y", "Vy_y", "Vz_y", "index_y",
                 1        9       0       1018        112       1        9       0       1018
                 2        8       0       1019        112       1        9       0       1018
                 1        9       0       1018        112       2        8       0       1019
                 2        8       0       1019        112       2        9       0       1019
        This dataframe represent ALL the possible atom couple (twice). However, since the operator value we want to compute is a+a+aa, we must delete rows where the left atom is the same than the right atom.
        Step 2 consists in computing simply the difference in each axis [Vx, Vy, Vz] and suppress all other columns. With our exemple, we would end up with
        atXY = "Vx"  "Vy"  "Vz"
                -1     1     0
                 1    -1     0
        We finally compute the 3D histogram of this matrix, using the properties of the voxel map entered by the user.

        Parameters
        ----------
        atX : pd.DataFrame
            Pandas DataFrame with atom speeds : it must contains "Vx", "Vy" and "Vz" and "index" where Vi are the speed along each axis and index is a unique number per atom of the dataset. Such df are created in the update_atoms_in_beams() method.
        atY : pd.DataFrame
            Pandas DataFrame with atoms of a given.

        Returns
        -------
        pd.DataFrame
            Pandas DataFrame with all velocities differences.
        """
        r = [
            self.voxel_map_range[0][0],
            self.voxel_map_range[0][1],
            self.voxel_map_range[1][0],
            self.voxel_map_range[1][1],
            self.voxel_map_range[2][0],
            self.voxel_map_range[2][1],
        ]
        # STEP 1 : create a "big" dataframe with all (ki,kj) couples with ki different from kj, ki is taken from atX and kj is taken from atY.
        if numerator is True:
            atXY = atX.merge(atY, how="outer", on="Cycle")
            # atXY = atX.merge(atY, how="cross") --> not possible with variable cycles
            atXY.drop(
                atXY[atXY["index_x"] == atXY["index_y"]].index,
                axis=0,
                inplace=True,
            )  # suppress identical atoms because a+a+aa is zero when acting on the SAME atom.
            # STEP 2 : compute difference or sum of the two atoms
            for axis in self.axis:
                if local is True:
                    atXY[axis] = atXY[axis + "_x"] - atXY[axis + "_y"]
                else:
                    atXY[axis] = atXY[axis + "_x"] + atXY[axis + "_y"]
            # STEP 3 : 3D histogram using torch.
            # atXY = atXY.astype(np.float32)
            G2XY, edges = torch.histogramdd(
                torch.from_numpy(np.array([atXY[ax].to_numpy() for ax in self.axis]).T),
                bins=list(self.voxel_map_size),
                range=r,
            )
            del atXY
            return G2XY.detach().cpu().numpy()
        elif numerator is False:
            ### === POSSIBILITY 1 : USE TORCH ===
            # STEP 1: create a really really big dataframe of size (big, 3), initialize as empty
            atXY = torch.zeros((len(atX) * len(atY), len(self.axis)))
            # STEP 2 : compute difference or sum and add it to the torch
            for i, ax in enumerate(self.axis):
                a = torch.cartesian_prod(
                    torch.tensor(atX[ax].to_numpy()), torch.tensor(atY[ax].to_numpy())
                )
                if local is True:
                    atXY[:, i] = a[:, 0] - a[:, 1]
                else:
                    atXY[:, i] = a[:, 0] + a[:, 1]
            # STEP 3 : 3D histogram using torch.
            G2XY, edges = torch.histogramdd(
                atXY, bins=list(self.voxel_map_size), range=r
            )
            ### === POSSIBILITY 2 : USE PANDAS ===
            # atXY = atX.merge(atY, how="cross")
            # atXY.drop(
            #     atXY[atXY["index_x"] == atXY["index_y"]].index,
            #     axis=0,
            #     inplace=True,
            # )  # suppress identical atoms because a+a+aa is zero when acting on the SAME atom.
            # # STEP 2 : compute difference or sum of the two atoms
            # for axis in self.axis:
            #     if local is True:
            #         atXY[axis] = atXY[axis + "_x"] - atXY[axis + "_y"]
            #     else:
            #         atXY[axis] = atXY[axis + "_x"] + atXY[axis + "_y"]
            return G2XY.detach().cpu().numpy()

    def compute_numerator(self):
        cycles_array_splitted = np.array_split(
            self.cycles_array, int(self.n_cycles / self.computer_performance)
        )
        for cycles in cycles_array_splitted:
            ### Beam A
            atA = self.atomsA[self.atomsA["Cycle"].isin(cycles)]
            G2AA = self.get_G2(atA, atA, local=True, numerator=True)
            self.result["G2AA"] += G2AA.flatten()
            # Beam B
            if self.only_one_beam is False:
                atB = self.atomsB[self.atomsB["Cycle"].isin(cycles)]
                G2BB = self.get_G2(atB, atB, local=True, numerator=True)
                self.result["G2BB"] += G2BB.flatten()
                ### Crossed A & B
                G2AB = self.get_G2(atA, atB, local=False, numerator=True)
                self.result["G2AB"] += G2AB.flatten()
            else:
                G2AA_crossed = self.get_G2(atA, atA, local=False, numerator=True)
                self.result["G2AB"] += G2AA_crossed.flatten()

    def compute_denominator(self):
        for cycle in tqdm(self.cycles_array):
            ### Beam A
            atA = self.atomsA[self.atomsA["Cycle"] == cycle]
            atA_all = self.atomsA[~(self.atomsA["Cycle"] == cycle)].sample(
                frac=self.denominator_ratio, replace=False
            )

            G2AA1 = self.get_G2(atA, atA_all, local=True, numerator=False)
            G2AA2 = self.get_G2(atA_all, atA, local=True, numerator=False)
            self.result["G2AA denominator"] += (G2AA1.flatten() + G2AA2.flatten()) / 2
            ### Beam B
            if self.only_one_beam is False:
                atB = self.atomsB[self.atomsB["Cycle"] == cycle]
                atB_all = self.atomsB[~(self.atomsB["Cycle"] == cycle)].sample(
                    frac=self.denominator_ratio, replace=False
                )
                G2BB1 = self.get_G2(atB, atB_all, local=True, numerator=False)
                G2BB2 = self.get_G2(atB_all, atB, local=True, numerator=False)
                self.result["G2BB denominator"] += (
                    G2BB1.flatten() + G2BB2.flatten()
                ) / 2
            ### Crossed Beam
            if self.only_one_beam is False:
                G2AB = self.get_G2(atA, atB_all, local=False, numerator=False)
                G2BA = self.get_G2(atB, atA_all, local=False, numerator=False)
                self.result["G2AB denominator"] += (G2AB.flatten() + G2BA.flatten()) / 2
            else:
                G2AB = self.get_G2(atA, atA_all, local=False, numerator=False)
            # self.result["G2AB denominator1"] += G2AB.flatten()
            # self.result["G2AB denominator2"] += G2BB.flatten()
        for G2 in ["G2AA", "G2BB", "G2AB"]:
            self.result[G2 + " denominator"] = self.result[G2 + " denominator"] / (
                self.n_cycles - 1
            )
        # self.result["G2AB denominator"]

    def compute_correlations(
        self,
    ):
        """This function initialize"""
        self.initialize_voxel_map_properties()
        self.update_atoms_in_beams()
        self.compute_numerator()
        self.compute_denominator()

    def compute_correlations_bootstrap(self, N=30):
        if "result" not in self.__dict__:
            self.initialize_voxel_map_properties()
            self.update_atoms_in_beams()
            self.compute_numerator()
            self.compute_denominator()

        if "G2AA mean" not in self.result.columns:
            self.boostrap_iteration = 0
            for X in ["G2AA", "G2BB", "G2AB"]:
                self.result[X + " mean"] = np.zeros(len(self.result))
                self.result[X + " squared"] = np.zeros(len(self.result))
                self.result[X + " std"] = np.zeros(len(self.result))

        ## Start to Bootstrap Atoms
        self.save_copy_of_atoms()
        for boostrap_num in tqdm(range(N)):
            self.bootstrap_atoms()
            self.update_atoms_in_beams()
            cycles_array_splitted = np.array_split(
                self.cycles_array, int(self.n_cycles / self.computer_performance)
            )
            res = pd.DataFrame(
                data=np.zeros((len(self.result), 3)), columns=["G2AA", "G2BB", "G2AB"]
            )

            for cycles in cycles_array_splitted:
                ### Beam A
                atA = self.atomsA[self.atomsA["Cycle"].isin(cycles)]
                G2AA = self.get_G2(atA, atA, local=True, numerator=True)
                res["G2AA"] += G2AA.flatten()
                atB = self.atomsB[self.atomsB["Cycle"].isin(cycles)]
                G2BB = self.get_G2(atB, atB, local=True, numerator=True)
                res["G2BB"] += G2BB.flatten()
                ### Crossed A & B
                G2AB = self.get_G2(atA, atB, local=False, numerator=True)
                res["G2AB"] += G2AB.flatten()
            # Save the result in result
            for X in ["G2AA", "G2BB", "G2AB"]:
                self.result[X + " mean"] = (
                    self.result[X + " mean"] * self.boostrap_iteration + res[X]
                ) / (self.boostrap_iteration + 1)
                self.result[X + " squared"] = (
                    self.result[X + " squared"] * self.boostrap_iteration + res[X] ** 2
                ) / (self.boostrap_iteration + 1)
            self.boostrap_iteration += 1
            for X in ["G2AA", "G2BB", "G2AB"]:
                self.result[X + " std"] = (
                    self.result[X + " squared"] - self.result[X + " mean"] ** 2
                )

    ################################################
    ## VISUALISATION FUNCTIONS CALLED BY THE USER ##
    ################################################
    def get_memory_informations(self):
        print("#" * 50)
        print("==== Default memory usage with initial data ====")
        mem = self.atoms.memory_usage(deep=True).sum() / 1024**2
        average_cycle_memory = mem / self.n_cycles
        a = self.atoms.groupby("Cycle").count()
        average_atom_number = a.mean()["Vx"]
        max_atom_number = np.max(a.mean()["Vx"])
        no_of_atoms = len(self.atoms)
        cycle_memory = mem / no_of_atoms * average_atom_number
        print("Memory needed for all atoms : {:.3f} MB".format(mem))
        print(
            "Memory needed to compute the denominator for the average cycle :  {:.3f} MB".format(
                mem * average_atom_number
            )
        )
        print(
            "Memory needed to compute the denominator for the worst cycle :  {:.3f} MB".format(
                mem * max_atom_number
            )
        )
        print(
            "Memory needed to compute the numerator with your computer_performance of {} : {:.3f} MB".format(
                self.computer_performance,
                cycle_memory * average_atom_number * self.computer_performance,
            )
        )

        self.update_atoms_in_beams()
        df = self.atomsA
        for col in self.atomsA.columns:
            print("Type of axis {} : {}".format(col, df[col].dtype))
        print("==== Memory usage with beam A ====")

        mem = df.memory_usage(deep=True).sum() / 1024**2
        average_cycle_memory = mem / self.n_cycles
        a = df.groupby("Cycle").count()
        average_atom_number = a.mean()["Vx"]
        max_atom_number = np.max(a.mean()["Vx"])
        no_of_atoms = len(self.atoms)
        cycle_memory = mem / no_of_atoms * average_atom_number
        print("Memory needed for all atoms : {:.3f} MB".format(mem))
        print(
            "Memory needed to compute the denominator for the average cycle :  {:.3f} MB".format(
                mem * average_atom_number
            )
        )
        print(
            "Memory needed to compute the denominator for the worst cycle :  {:.3f} MB".format(
                mem * max_atom_number
            )
        )
        print(
            "Memory needed to compute the numerator with your computer_performance of {} : {:.3f} MB".format(
                self.computer_performance,
                cycle_memory * average_atom_number * self.computer_performance,
            )
        )

    def show_density_plots(self):
        # First Plot : 2D density
        fig, axes = plt.subplots(figsize=(12, 4), ncols=3)
        for i in range(3):
            x = self.axis[i]
            y = self.axis[(i + 1) % 3]
            ax = axes[i]
            sns.histplot(
                self.atoms, x=x, y=y, ax=axes[i], cbar=False, cmap=plt.cm.coolwarm
            )  # , palette = "twilight")

            def draw_box(cX, σX, cY, σY, **kwargs):
                ax.plot(
                    [cX - σX, cX + σX, cX + σX, cX - σX, cX - σX],
                    [cY - σY, cY - σY, cY + σY, cY + σY, cY - σY],
                    **kwargs,
                )

            # on affiche la boite 1
            posX = self.beams["A"][x]["position"]
            sizeX = self.beams["A"][x]["size"]
            posY = self.beams["A"][y]["position"]
            sizeY = self.beams["A"][y]["size"]
            draw_box(
                posX, sizeX / 2, posY, sizeY / 2, color="darkgreen", label="beam A"
            )
            # draw_box(posX, 3*sizeX / 2, posY, 3*sizeY / 2, color="darkgreen",ls = "--")
            # on affiche la boite 1
            posX = self.beams["B"][x]["position"]
            sizeX = self.beams["B"][x]["size"]
            posY = self.beams["B"][y]["position"]
            sizeY = self.beams["B"][y]["size"]
            draw_box(posX, sizeX / 2, posY, sizeY / 2, color="dimgray", label="beam B")
            # draw_box(posX, 3*sizeX / 2, posY, 3*sizeY / 2, color="dimgray",ls = "--")
        plt.legend(loc=1)
        plt.tight_layout()
        plt.show()

        # Second Plot : every thing in ROI
        fig, axes = plt.subplots(figsize=(12, 4), ncols=3)
        for i in range(3):
            x = self.axis[i]
            ax = axes[i]
            sns.histplot(data=self.atoms, x=x, ax=ax, color="C" + str(i))
        plt.tight_layout()
        plt.show()
        # third plot : each beam
        self.update_atoms_in_beams()

        fig, axes = plt.subplots(figsize=(12, 4), ncols=3)
        statsA = []
        statsB = []
        for i in range(3):
            # sns.histplot(data=self.atoms, x =self.axis[i] , ax = axes[i],
            #               color = "C2", alpha = 0.2, label = "ROI")
            sns.histplot(
                data=self.atomsA,
                x=self.axis[i],
                ax=axes[i],
                color="C0",
                alpha=0.5,
                label="beam A",
            )
            sns.histplot(
                data=self.atomsB,
                x=self.axis[i],
                ax=axes[i],
                color="C1",
                alpha=0.5,
                label="beam B",
                kde=True,
            )
            statsA.append(str(self.atomsA[self.axis[i]].describe()))
        plt.legend()
        plt.tight_layout()
        plt.show()
        print("\n " + "=" * 30 + " BEAM A " + "=" * 30 + "\n")
        for elem, val in self.beams["A"].items():
            print(elem, val)
        print(self.atomsA.describe())
        print("\n " + "=" * 30 + " BEAM B " + "=" * 30 + "\n")
        for elem, val in self.beams["B"].items():
            print(elem, val)
        print(self.atomsB.describe())


class CorrelationHe2StyleBigDenominator(CorrelationHe2Style):
    """
    Classe selfelationHe2Style. Prend en entrée un dataframe avec X, Y et T et calcule de corrélations etc...
    """

    def __init__(self, atoms, **kwargs):
        """
        Object initialization, sets parameters as the user defined, build the atoms dataframe and apply ROD and ROI.
        """
        super().__init__(atoms, **kwargs)

    def get_G2(
        self, atX: pd.DataFrame, atY: pd.DataFrame, local=True, numerator=True
    ) -> pd.DataFrame:
        """Function that compute the 3D velocity difference between all atoms in the crossed atX x atY dataframe (atoms in beam X, X being A or B). If local is True, it computes the difference while if local is False, it computes the sum.

        Detailed explanation :
        atX = DataFrame("Vx":[1,2] ; "Vy":[9,8] ; "Vz":[0,0], "index": [1018, 1019], "Cycle":[112,112] )
        atY = atX and we want to compute the local correlation. First step is to create the cross product of this DataFrame using merge (or concat depending on the version of the code). We create the following dataframe (_x and_y suffix are added by pandas) :
        atXY = "Vx_x", "Vy_x", "Vz_x", "index_x", "Cycle_x", "Vx_y", "Vy_y", "Vz_y", "index_y", "Cycle_y")
                 1        9       0       1018        112       1        9       0       1018        112
                 2        8       0       1019        112       1        9       0       1018        112
                 1        9       0       1018        112       2        8       0       1019        112
                 2        8       0       1019        112       2        9       0       1019        112
        This dataframe represent ALL the possible atom couple (twice). However, since the operator value we want to compute is a+a+aa, we must delete rows where the left atom is the same than the right atom.
        We then compute simply the difference in each axis [Vx, Vy, Vz] and suppress all other columns. With our exemple, we would end up with
        atXY = "Vx"  "Vy"  "Vz"
                -1     1     0
                 1    -1     0
        We finally compute the 3D histogram of this matrix, using the properties of the voxel map entered by the user.

        Parameters
        ----------
        atX : pd.DataFrame
            Pandas DataFrame with atom speeds : it must contains "Vx", "Vy" and "Vz" and "index" where Vi are the speed along each axis and index is a unique number per atom of the dataset. Such df are created in the update_atoms_in_beams() method.
        atY : pd.DataFrame
            Pandas DataFrame with atoms of a given.

        Returns
        -------
        pd.DataFrame
            Pandas DataFrame with all velocities differences.
        """
        if numerator:
            atXY = atX.merge(atY, how="outer", on="Cycle")
            atXY.drop(
                atXY[atXY["index_x"] == atXY["index_y"]].index,
                axis=0,
                inplace=True,
            )
        else:
            atXY = atX.merge(atY, how="cross")
            # suppress identical atoms because a+a+aa is zero when acting on the SAME atom.
        for axis in self.axis:
            if local:
                atXY[axis] = atXY[axis + "_x"] - atXY[axis + "_y"]
            else:
                if axis == "Vz":
                    atXY[axis] = atXY[axis + "_x"] + atXY[axis + "_y"]
                else:
                    atXY[axis] = atXY[axis + "_x"] + atXY[axis + "_y"]

        cols_to_remove = [col for col in atXY.columns if col not in self.axis]
        # print("Datasize : {}".format(len(atXY)))
        r = [
            self.voxel_map_range[0][0],
            self.voxel_map_range[0][1],
            self.voxel_map_range[1][0],
            self.voxel_map_range[1][1],
            self.voxel_map_range[2][0],
            self.voxel_map_range[2][1],
        ]
        atXY.drop(cols_to_remove, axis=1, inplace=True)
        # t0 = time.time()
        # G2XY, edges = np.histogramdd(
        #    atXY.to_numpy(), bins=self.voxel_map_size, range=self.voxel_map_range
        # )
        # t1 = time.time()
        G2XY, edges = torch.histogramdd(
            torch.from_numpy(atXY.to_numpy()), bins=list(self.voxel_map_size), range=r
        )
        # t2 = time.time()
        # dt_torch = t2 - t1
        # dt_numpy = t1 - t0
        # if dt_torch < dt_numpy:
        #     print("Torch is {:.0f} times better".format((dt_numpy/dt_torch - 1)))
        # else:
        #     print("Numpy is {:.0f} times better".format((dt_torch/dt_numpy - 1)))
        # del atXY
        return G2XY.detach().cpu().numpy()

    def get_all_atoms_except_cycle(self, atoms, cycles):
        return atoms[~(atoms["Cycle"].isin(cycles))]

    def compute_correlations(self):
        """Principal function of the code. L'idée du code est le suivant :
        * on parcourt les cycles de la séquence et pour chaque cycle, on récupère les atomes dans la ROI1 et la ROI2.
        * pour chaque cycle, on calcule la différence (corrélations locales) ou la somme (corrélation croisées) de l'ensemble des couple d'atomes.
        * on fait ensuite l'histogramme 3D de cet ensemble. On l'ajoute à l'histogramme total.
        """
        self.initialize_voxel_map_properties()
        self.update_atoms_in_beams()

        cycles_array_splitted = np.array_split(
            self.cycles_array, int(self.n_cycles / self.computer_capacity)
        )
        self.normalization_factors = []
        for i, cycles in tqdm(enumerate(cycles_array_splitted)):
            normalization_factor = self.n_cycles + 1 - len(cycles)
            self.normalization_factors.append(normalization_factor)
            ### Beam A
            atA = self.atomsA[self.atomsA["Cycle"].isin(cycles)]
            atA_rand = self.get_all_atoms_except_cycle(self.atomsA, cycles)
            G2AA = self.get_G2(atA, atA, local=True, numerator=True)
            self.result["G2AA"] += G2AA.flatten()
            G2AA_rand = self.get_G2(atA, atA_rand, local=True, numerator=False)
            self.result["G2AA random"] += G2AA_rand.flatten() / normalization_factor
            if self.only_one_beam is False:
                ### Beam B
                atB = self.atomsB[self.atomsB["Cycle"].isin(cycles)]
                atB_rand = self.get_all_atoms_except_cycle(self.atomsB, cycles)
                G2BB = self.get_G2(atB, atB, local=True, numerator=True)
                self.result["G2BB"] += G2BB.flatten()
                G2BB_rand = self.get_G2(atB, atB_rand, local=True, numerator=False)
                self.result["G2BB random"] += G2BB_rand.flatten() / normalization_factor
                ### Crossed A & B
                G2AB = self.get_G2(atA, atB, local=False, numerator=True)
                self.result["G2AB"] += G2AB.flatten()
                G2AB_rand = self.get_G2(atA, atB_rand, local=False, numerator=False)
                self.result["G2AB random"] += G2AB_rand.flatten() / normalization_factor
            # if i ==3:
            #     break
        print(
            "[Info] : The denominator normalization used {:.1f} % of atoms".format(
                100 * np.mean(self.normalization_factors) / self.n_cycles
            )
        )
        self.result["g2 aa"] = self.result["G2AA"] / self.result["G2AA random"]
        # self.result["G2AA"] = self.beamA_volume * self.result["G2AA"] --> normalisation ? not needed yet but keep in mind
        # self.result["G2AA random"] = self.beamA_volume * self.result["G2AA"]
        if self.only_one_beam is False:
            self.result["g2 bb"] = self.result["G2BB"] / self.result["G2BB random"]
            self.result["g2 ab"] = self.result["G2AB"] / self.result["G2AB random"]


import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from random import choice, shuffle, randint, sample
from tqdm import tqdm, trange
import copy
import time
import torch
from math import ceil

import time


class CorrelationHe2StyleBigDenominatorTest:
    """
    Classe selfelationHe2Style. Prend en entrée un dataframe avec X, Y et T et calcule de corrélations etc...
    """

    def __init__(self, atoms, **kwargs):
        """
        Object initialization, sets parameters as the user defined, build the atoms dataframe and apply ROD and ROI.
        """
        # super().__init__(atoms, **kwargs)
        self.atoms = copy.deepcopy(atoms)
        self.ROI = {}
        self.ROD = {}
        self.round_decimal = 7
        self.only_one_beam = False
        self.computer_performance = 1
        self.computer_capacity = 1
        self.ratio_of_atoms_for_denominator = 0.1
        self.axis = ["Vx", "Vy", "Vz"]
        self.ref_frame_speed = {"Vx": 0, "Vy": 0, "Vz": 0}
        self.ROI = {}
        self.voxel_size = {"Vx": 0.1, "Vy": 0.1, "Vz": 0.1}
        self.voxel_numbers = {"Vx": 11, "Vy": 11, "Vz": 11}
        self.beams = {
            "A": {
                "Vx": {"size": 20, "position": 0},
                "Vy": {"size": 20, "position": 0},
                "Vz": {"size": 5, "position": 25},
            },
            "B": {
                "Vx": {"size": 20, "position": 0},
                "Vy": {"size": 20, "position": 0},
                "Vz": {"size": 5, "position": -25},
            },
        }
        self.__dict__.update(kwargs)
        # atoms : (mm,mm, ms) --> (mm/s, mm/s, mm/s)
        # self.build_the_atoms_dataframe()
        self.apply_ROI()  # Keep only atoms in ROI
        self.apply_ROD()  # Take off atoms in Region Of Desinterest.
        self.cycles_array = atoms["Cycle"].unique()
        self.n_cycles = len(atoms["Cycle"].unique())
        self.random_cycle_mapping = self.get_randomize_cycle_mapping()

    def get_G2(
        self, atX: pd.DataFrame, atY: pd.DataFrame, local=True, numerator=True
    ) -> pd.DataFrame:
        """Function that compute the 3D velocity difference between all atoms in the crossed atX x atY dataframe (atoms in beam X, X being A or B). If local is True, it computes the difference while if local is False, it computes the sum.

        Detailed explanation :
        atX = DataFrame("Vx":[1,2] ; "Vy":[9,8] ; "Vz":[0,0], "index": [1018, 1019], "Cycle":[112,112] )
        atY = atX and we want to compute the local correlation. First step is to create the cross product of this DataFrame using merge (or concat depending on the version of the code). We create the following dataframe (_x and_y suffix are added by pandas) :
        atXY = "Vx_x", "Vy_x", "Vz_x", "index_x", "Cycle_x", "Vx_y", "Vy_y", "Vz_y", "index_y", "Cycle_y")
                 1        9       0       1018        112       1        9       0       1018        112
                 2        8       0       1019        112       1        9       0       1018        112
                 1        9       0       1018        112       2        8       0       1019        112
                 2        8       0       1019        112       2        9       0       1019        112
        This dataframe represent ALL the possible atom couple (twice). However, since the operator value we want to compute is a+a+aa, we must delete rows where the left atom is the same than the right atom.
        We then compute simply the difference in each axis [Vx, Vy, Vz] and suppress all other columns. With our exemple, we would end up with
        atXY = "Vx"  "Vy"  "Vz"
                -1     1     0
                 1    -1     0
        We finally compute the 3D histogram of this matrix, using the properties of the voxel map entered by the user.

        Parameters
        ----------
        atX : pd.DataFrame
            Pandas DataFrame with atom speeds : it must contains "Vx", "Vy" and "Vz" and "index" where Vi are the speed along each axis and index is a unique number per atom of the dataset. Such df are created in the update_atoms_in_beams() method.
        atY : pd.DataFrame
            Pandas DataFrame with atoms of a given.

        Returns
        -------
        pd.DataFrame
            Pandas DataFrame with all velocities differences.
        """
        if numerator:
            atXY = atX.merge(atY, how="outer", on="Cycle")
            atXY.drop(
                atXY[atXY["index_x"] == atXY["index_y"]].index,
                axis=0,
                inplace=True,
            )
        else:
            atXY = atX.merge(atY, how="cross")
            # suppress identical atoms because a+a+aa is zero when acting on the SAME atom.
        for axis in self.axis:
            if local:
                atXY[axis] = atXY[axis + "_x"] - atXY[axis + "_y"]
            else:
                atXY[axis] = atXY[axis + "_x"] + atXY[axis + "_y"]
        cols_to_remove = [col for col in atXY.columns if col not in self.axis]
        # print("Datasize : {}".format(len(atXY)))
        r = [
            self.voxel_map_range[0][0],
            self.voxel_map_range[0][1],
            self.voxel_map_range[1][0],
            self.voxel_map_range[1][1],
            self.voxel_map_range[2][0],
            self.voxel_map_range[2][1],
        ]
        atXY.drop(cols_to_remove, axis=1, inplace=True)
        # t0 = time.time()
        # G2XY, edges = np.histogramdd(
        #    atXY.to_numpy(), bins=self.voxel_map_size, range=self.voxel_map_range
        # )
        # t1 = time.time()
        G2XY, edges = torch.histogramdd(
            torch.from_numpy(atXY.to_numpy()), bins=list(self.voxel_map_size), range=r
        )
        # t2 = time.time()
        # dt_torch = t2 - t1
        # dt_numpy = t1 - t0
        # if dt_torch < dt_numpy:
        #     print("Torch is {:.0f} times better".format((dt_numpy/dt_torch - 1)))
        # else:
        #     print("Numpy is {:.0f} times better".format((dt_torch/dt_numpy - 1)))
        # del atXY
        return G2XY.detach().cpu().numpy()

    def get_all_atoms_except_cycle(self, atoms, cycles):
        at = atoms[~(atoms["Cycle"].isin(cycles))]
        if self.ratio_of_atoms_for_denominator == 1:
            return at, self.n_cycles + 1 - len(cycles)
        possible_cycles = at["Cycle"].unique()
        no_of_cycles = ceil(len(possible_cycles) * self.ratio_of_atoms_for_denominator)
        selected_cycles = sample(list(possible_cycles), no_of_cycles)
        at = at[(at["Cycle"].isin(selected_cycles))]
        return at, no_of_cycles

    def compute_correlations(self):
        """Principal function of the code. L'idée du code est le suivant :
        * on parcourt les cycles de la séquence et pour chaque cycle, on récupère les atomes dans la ROI1 et la ROI2.
        * pour chaque cycle, on calcule la différence (corrélations locales) ou la somme (corrélation croisées) de l'ensemble des couple d'atomes.
        * on fait ensuite l'histogramme 3D de cet ensemble. On l'ajoute à l'histogramme total.
        """
        self.initialize_voxel_map_properties()
        self.update_atoms_in_beams()

        cycles_array_splitted = np.array_split(
            self.cycles_array, int(self.n_cycles / self.computer_capacity)
        )
        self.normalization_factors = []
        for i, cycles in enumerate(cycles_array_splitted):
            # normalization_factor = self.n_cycles + 1 - len(cycles)
            # self.normalization_factors.append(normalization_factor)
            ### Beam A
            self.time = time.time()
            atA = self.atomsA[self.atomsA["Cycle"].isin(cycles)]
            atA_rand, normalization_factor = self.get_all_atoms_except_cycle(
                self.atomsA, cycles
            )
            print("Gather data for one cycle : " + str(int(time.time() - self.time)))
            self.time = time.time()
            G2AA = self.get_G2(atA, atA, local=True, numerator=True)
            print("Compute G2 numerator : " + str(int(time.time() - self.time)))
            self.time = time.time()
            self.result["G2AA"] += G2AA.flatten()
            G2AA_rand = self.get_G2(atA, atA_rand, local=True, numerator=False)
            print("Compute G2 denominator : " + str(int(time.time() - self.time)))
            self.time = time.time()
            self.result["G2AA random"] += G2AA_rand.flatten() / normalization_factor
            if self.only_one_beam is False:
                ### Beam B
                atB = self.atomsB[self.atomsB["Cycle"].isin(cycles)]
                atB_rand = self.get_all_atoms_except_cycle(self.atomsB, cycles)
                G2BB = self.get_G2(atB, atB, local=True, numerator=True)
                self.result["G2BB"] += G2BB.flatten()
                G2BB_rand = self.get_G2(atB, atB_rand, local=True, numerator=False)
                self.result["G2BB random"] += G2BB_rand.flatten() / normalization_factor
                ### Crossed A & B
                G2AB = self.get_G2(atA, atB, local=False, numerator=True)
                self.result["G2AB"] += G2AB.flatten()
                G2AB_rand = self.get_G2(atA, atB_rand, local=False, numerator=False)
                self.result["G2AB random"] += G2AB_rand.flatten() / normalization_factor
            # if i ==3:
            #     break
        print(
            "[Info] : The denominator normalization used {:.1f} % of atoms".format(
                100 * np.mean(self.normalization_factors) / self.n_cycles
            )
        )
        self.result["g2 aa"] = self.result["G2AA"] / self.result["G2AA random"]
        # self.result["G2AA"] = self.beamA_volume * self.result["G2AA"] --> normalisation ? not needed yet but keep in mind
        # self.result["G2AA random"] = self.beamA_volume * self.result["G2AA"]
        if self.only_one_beam is False:
            self.result["g2 bb"] = self.result["G2BB"] / self.result["G2BB random"]
            self.result["g2 ab"] = self.result["G2AB"] / self.result["G2AB random"]

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

    def initialize_voxel_map_properties(self):
        """_summary_"""
        self.voxel_map_size = tuple(
            [int(self.voxel_numbers[axis]) for axis in self.axis]
        )
        # self.voxel_map = np.zeros(self.voxel_map_size)
        # self.voxel_size = {"Vx": 0.1, "Vy": 0.1, "Vz": 0.1}
        self.voxel_map_range = tuple(
            [
                (
                    -self.voxel_numbers[axis] * self.voxel_size[axis] / 2,
                    self.voxel_numbers[axis] * self.voxel_size[axis] / 2,
                )
                for axis in self.axis
            ]
        )
        voxel_centers = []
        for axis in self.axis:
            mini = (
                -self.voxel_numbers[axis] * self.voxel_size[axis] / 2
                + self.voxel_size[axis] / 2
            )
            maxi = +self.voxel_numbers[axis] * self.voxel_size[axis] / 2
            voxel_centers.append(np.arange(mini, maxi, step=self.voxel_size[axis]))
            if len(voxel_centers[-1]) != self.voxel_numbers[axis]:
                print(len(voxel_centers[-1]))
                raise Exception(
                    "Something strange appended in the voxel map initialization."
                )
        # create the result dataframe

        data = np.array(
            [
                [x, y, z]
                for x in voxel_centers[0]
                for y in voxel_centers[1]
                for z in voxel_centers[2]
            ]
        )
        self.result = pd.DataFrame(data=data, columns=self.axis)
        if self.only_one_beam is True:
            for column in [
                "G2AA",
                "G2AA random",
            ]:
                self.result[column] = np.zeros(len(self.result))
        else:
            for column in [
                "G2AA",
                "G2BB",
                "G2AB",
                "G2AA random",
                "G2BB random",
                "G2AB random1",
                "G2AB random2",
                "G2AB random",
            ]:
                self.result[column] = np.zeros(len(self.result))

    def get_atoms_in_beams(self, beam):
        """
        Récupère  les atomes à l'intérieur de beam. Cela permet d'alléger les données à traiter.
        Si le beam est vide, la méthode ne fait rien. Le format du beam doit être {"Vx": {"position":120, "size":-120}}
        """
        if beam:
            atoms_in_beam = self.atoms
            for key, entry in beam.items():
                atoms_in_beam = atoms_in_beam[
                    (
                        (atoms_in_beam[key] <= entry["position"] + entry["size"] / 2.0)
                        & (atoms_in_beam[key] > entry["position"] - entry["size"] / 2.0)
                    )
                ]
        return atoms_in_beam

    def update_atoms_in_beams(self):
        self.atomsA = self.get_atoms_in_beams(beam=self.beams["A"]).reset_index()
        self.atomsA["index"] = np.arange(0, len(self.atomsA))
        self.atomsB = self.get_atoms_in_beams(beam=self.beams["B"]).reset_index()
        self.atomsB["index"] = np.arange(
            len(self.atomsA), len(self.atomsA) + len(self.atomsB)
        )
        self.atomsA_randomized = self.atomsA.copy()
        self.atomsA_randomized["Cycle"] = self.atomsA_randomized["Cycle"].map(
            self.random_cycle_mapping
        )
        self.atomsB_randomized = self.atomsB.copy()
        self.atomsB_randomized["Cycle"] = self.atomsB_randomized["Cycle"].map(
            self.random_cycle_mapping
        )
        self.beamA_volume = 1
        self.beamB_volume = 1
        for ax in self.axis:
            self.beamA_volume *= self.beams["A"][ax]["size"]
            self.beamB_volume *= self.beams["B"][ax]["size"]

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

    def get_randomize_cycle_mapping(self) -> dict:
        """This function shuffles the self.cycles_array to create a random_cycles array. It ensure that all cycles have been mooved and then defined a mapping between the two arrays.

        Returns
        -------
        dictionary
            mapping function to map initial cycle to randomized cycle
        """
        random_cycles = self.cycles_array.copy()
        shuffle(random_cycles)
        while np.min((random_cycles - self.cycles_array) ** 2) == 0:
            a = randint(0, self.n_cycles)
            arg = np.argmin((random_cycles - self.cycles_array) ** 2)
            random_cycles[arg], random_cycles[(arg + a) % self.n_cycles] = (
                random_cycles[(arg + a) % self.n_cycles],
                random_cycles[arg],
            )
        return dict(zip(self.cycles_array, random_cycles))


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
        "Vx": {"min": -90, "max": 90},
        "Vy": {"min": -90, "max": 90},
        "Vz": {"min": -50, "max": 50},
    }
    posiZ = 22
    REF_FRAME_SPEED = {"Vx": 0, "Vy": -4.3, "Vz": 94}
    ARRIVAL_TIME = 307.55450877426836
    beams = {
        "A": {
            "Vx": {"size": 100, "position": 0},
            "Vy": {"size": 100, "position": 0},
            "Vz": {"size": 25, "position": -posiZ},
        },
        "B": {
            "Vx": {"size": 100, "position": 0},
            "Vy": {"size": 100, "position": 0},
            "Vz": {"size": 25, "position": posiZ},
        },
    }
    voxel_size = {"Vx": 2, "Vy": 2, "Vz": 0.25}
    computer_capacity = 300
    # selected_data = selected_data[selected_data["Cycle"] < selected_data["Cycle"].unique()[300]]
    voxel_numbers = {"Vx": 61, "Vy": 61, "Vz": 151}
    selec_bec_arrival_times["BEC Arrival Time"] = selec_bec_arrival_times[
        "BEC Arrival Time with fit"
    ]
    corr = CorrelationHe2Style(
        selected_data,
        ROI=ROI,
        bec_arrival_time=ARRIVAL_TIME,
        ref_frame_speed=REF_FRAME_SPEED,
        beams=beams,
        voxel_size=voxel_size,
        voxel_numbers=voxel_numbers,
        computer_performance=computer_capacity,
    )
    # corr.atoms["Vz"] = corr.atoms["Vz"] - 0.55
    # corr.show_density_plots()
    corr.get_memory_informations()
    corr.compute_correlations()
    corr.result.to_pickle(
        "/home/victor/ownCloud/LabWiki/Journal/2023/06/13/data_corr.pkl"
    )
