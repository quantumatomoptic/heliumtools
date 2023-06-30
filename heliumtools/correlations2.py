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
import snoop
from random import choice, shuffle, randint
from tqdm import tqdm, trange
import copy
import time
import torch


class CorrelationHe2Style:
    """
    Classe selfelationHe2Style. Prend en entrée un dataframe avec X, Y et T et calcule de corrélations etc...
    """

    def __init__(self, atoms, **kwargs):
        """
        Object initialization, sets parameters as the user defined, build the atoms dataframe and apply ROD and ROI.
        """
        self.atoms = copy.deepcopy(atoms)
        self.bec_arrival_time = 307.763  # temps d'arrivée du BEC, en ms
        self.theoretical_arrival_time = 307.763  # 24/06/2022 & 17/05/2022
        self.raman_kick = 42.5  # mm/s, kick Raman
        self.gravity = 9.81
        self.bad_shot_limit = 100
        self.remove_shot_noise = True
        self.additionnal_bec_speed_ms = 0
        self.ROI = {}
        self.ROD = {}
        self.round_decimal = 7
        self.only_one_beam = False
        self.computer_performance = 1
        self.id = int(time.time())
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
        self.build_the_atoms_dataframe()
        self.apply_ROI()  # Keep only atoms in ROI
        self.apply_ROD()  # Take off atoms in Region Of Desinterest.
        self.cycles_array = atoms["Cycle"].unique()
        self.n_cycles = len(atoms["Cycle"].unique())
        self.random_cycle_mapping = self.get_randomize_cycle_mapping()
        # Initialisation de self.result datframe avec toutes les corrélations
        # self.result = pd.DataFrame(
        #     np.zeros(1, len(self.quantity_of_interest())),
        #     columns=self.quantity_of_interest,
        # )

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
            print(self.bec_arrival_time["BEC Arrival Time"].mean())
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
            print(self.bec_arrival_time["BEC Arrival Time"].mean())
            if "BEC Center X" in self.bec_arrival_time.columns:
                self.bec_arrival_time["BEC X"] = (
                    1000
                    * self.bec_arrival_time["BEC Center X"]
                    / self.bec_arrival_time["BEC Arrival Time"]
                )
            else:
                self.bec_arrival_time["BEC X"] = self.raman_kick
            print(self.bec_arrival_time["BEC Arrival Time"].mean())
            if "BEC Center Y" in self.bec_arrival_time.columns:
                self.bec_arrival_time["BEC Y"] = (
                    1000
                    * self.bec_arrival_time["BEC Center Y"]
                    / self.bec_arrival_time["BEC Arrival Time"]
                )
            else:
                self.bec_arrival_time["BEC Y"] = 0
            print(self.bec_arrival_time["BEC Arrival Time"].mean())
            self.bec_arrival_time["BEC Z"] = (
                0.5 * self.gravity * self.bec_arrival_time["BEC Arrival Time"]
                - l_fall / self.bec_arrival_time["BEC Arrival Time"]
            )
            # check now if we have perform some Bragg diffraction to fit the BEC
            print(self.bec_arrival_time["BEC Arrival Time"].mean())

            self.atoms = self.merge_dataframe_on_cycles(
                self.atoms, self.bec_arrival_time
            )
            # take off the speed of each BEC
            self.atoms["Vz"] = self.atoms["Vz"] - self.atoms["BEC Z"]
            self.atoms["Vx"] = self.atoms["Vx"] - self.atoms["BEC X"]
            self.atoms["Vy"] = self.atoms["Vy"] - self.atoms["BEC Y"]
            print(self.bec_arrival_time["BEC Arrival Time"].mean())

            ## On supprime ensuite les bad shots : là où il y a moins de 100 atomes
            df = self.atoms[self.atoms["Number of Atoms"] < self.bad_shot_limit]
            for cycle, nb_atoms in zip(df["Cycle"], df["Number of Atoms"]):
                print(f"Delete cycle {cycle} with only {nb_atoms} atoms.")
                self.n_cycles -= 1
            # drop columns with no more interest (for later calculations)
            for column in self.atoms.columns:
                if column not in ["Vz", "Vx", "Vy", "Cycle"]:
                    self.atoms.drop(column, inplace=True, axis=1)
            print(self.bec_arrival_time["BEC Arrival Time"].mean())

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

    def update_atoms_in_beams(self):
        self.atomsA = self.get_atoms_in_beams(beam=self.beams["A"]).reset_index()
        self.atomsA["index"] = np.arange(0, len(self.atomsA))
        self.atomsB = self.get_atoms_in_beams(beam=self.beams["B"]).reset_index()
        self.atomsB["index"] = np.arange(len(self.atomsA), len(self.atomsA)+len(self.atomsB))
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

    def get_G2(self, atX: pd.DataFrame, atY: pd.DataFrame, local=True) -> pd.DataFrame:
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
        atXY = atX.merge(
            atY, how="outer", on="Cycle"
        )
        # atXY = atX.merge(atY, how="cross")
        atXY.drop(
            atXY[atXY["index_x"] == atXY["index_y"]].index,
            axis=0,
            inplace=True,
        )  # suppress identical atoms because a+a+aa is zero when acting on the SAME atom.
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
        for cycles in tqdm(cycles_array_splitted):
            ### Beam A
            atA = self.atomsA[self.atomsA["Cycle"].isin(cycles)]
            atA_rand = self.atomsA_randomized[
                self.atomsA_randomized["Cycle"].isin(cycles)
            ]
            G2AA = self.get_G2(atA, atA, local=True)
            self.result["G2AA"] += G2AA.flatten()
            G2AA_rand = self.get_G2(atA, atA_rand, local=True)
            self.result["G2AA random"] += G2AA_rand.flatten()
            if self.only_one_beam is False:
                ### Beam B
                atB = self.atomsB[self.atomsB["Cycle"].isin(cycles)]
                atB_rand = self.atomsB_randomized[
                    self.atomsB_randomized["Cycle"].isin(cycles)
                ]
                G2BB = self.get_G2(atB, atB, local=True)
                self.result["G2BB"] += G2BB.flatten()
                G2BB_rand = self.get_G2(atB, atB_rand, local=True)
                self.result["G2BB random"] += G2BB_rand.flatten()
                ### Crossed A & B
                G2AB = self.get_G2(atA, atB, local=False)
                self.result["G2AB"] += G2AB.flatten()
                G2AB_rand1 = self.get_G2(atA, atB_rand, local=False)
                G2AB_rand2 = self.get_G2(atA_rand, atB, local=False)
                self.result["G2AB random1"] += G2AB_rand1.flatten()
                self.result["G2AB random2"] += G2AB_rand2.flatten()
                self.result["G2AB random2"] +=0.5*(G2AB_rand1.flatten() + G2AB_rand2.flatten())

        self.result["g2 aa"] = self.result["G2AA"] / self.result["G2AA random"]
        # self.result["G2AA"] = self.beamA_volume * self.result["G2AA"] --> normalisation ? not needed yet but keep in mind
        # self.result["G2AA random"] = self.beamA_volume * self.result["G2AA"]
        if self.only_one_beam is False:
            self.result["g2 bb"] = self.result["G2BB"] / self.result["G2BB random"]
            self.result["g2 ab"] = (
                2
                * self.result["G2AB"]
                / (self.result["G2AB random1"] + self.result["G2AB random2"])
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
        

    def get_G2(self, atX: pd.DataFrame, atY: pd.DataFrame, local=True, numerator = True) -> pd.DataFrame:
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
            atXY = atX.merge(
                atY, how="outer", on="Cycle"
            )
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
            G2AA = self.get_G2(atA, atA, local=True, numerator = True)
            self.result["G2AA"] += G2AA.flatten()
            G2AA_rand = self.get_G2(atA, atA_rand, local=True, numerator = False)
            self.result["G2AA random"] += G2AA_rand.flatten()/normalization_factor
            if self.only_one_beam is False:
                ### Beam B
                atB = self.atomsB[self.atomsB["Cycle"].isin(cycles)]
                atB_rand = self.get_all_atoms_except_cycle(self.atomsB, cycles)
                G2BB = self.get_G2(atB, atB, local=True, numerator = True)
                self.result["G2BB"] += G2BB.flatten()
                G2BB_rand = self.get_G2(atB, atB_rand, local=True, numerator = False)
                self.result["G2BB random"] += G2BB_rand.flatten()/normalization_factor
                ### Crossed A & B
                G2AB = self.get_G2(atA, atB, local=False, numerator = True)
                self.result["G2AB"] += G2AB.flatten()
                G2AB_rand = self.get_G2(atA, atB_rand, local=False, numerator = False)
                self.result["G2AB random"] += G2AB_rand.flatten()/normalization_factor
            # if i ==3:
            #     break
        print("[Info] : The denominator normalization used {:.1f} % of atoms".format(100 * np.mean(self.normalization_factors)/self.n_cycles))
        self.result["g2 aa"] = self.result["G2AA"] / self.result["G2AA random"]
        # self.result["G2AA"] = self.beamA_volume * self.result["G2AA"] --> normalisation ? not needed yet but keep in mind
        # self.result["G2AA random"] = self.beamA_volume * self.result["G2AA"]
        if self.only_one_beam is False:
            self.result["g2 bb"] = self.result["G2BB"] / self.result["G2BB random"]
            self.result["g2 ab"] = self.result["G2AB"] / self.result["G2AB random"]
            



if __name__ == "__main__":
    import seaborn as sns
    import pandas as pd
    selected_data = pd.read_pickle("/home/victor/ownCloud/LabWiki/Journal/2023/06/13/data.pkl")
    selec_bec_arrival_times = pd.read_pickle("/home/victor/ownCloud/LabWiki/Journal/2023/06/13/bec.pkl")
    ROI = {"Vx": {"min": -90, "max": 90},
        "Vy": {"min": -90, "max":90},
        "Vz": {"min": -50, "max": 50},
        }
    posiZ = 22
    REF_FRAME_SPEED = {"Vx": 0, "Vy": -4.3, "Vz": 94 }
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
    voxel_size = {"Vx": 2, "Vy": 2, "Vz": 0.1}
    computer_capacity = 8
    voxel_numbers = {"Vx": 31, "Vy":31 , "Vz": 101}
    selec_bec_arrival_times["BEC Arrival Time"] = selec_bec_arrival_times["BEC Arrival Time with fit"]
    corr = CorrelationHe2StyleBigDenominator(
        selected_data,
        ROI = ROI,
        bec_arrival_time=ARRIVAL_TIME,
        ref_frame_speed = REF_FRAME_SPEED,
        beams = beams, voxel_size = voxel_size,
        voxel_numbers = voxel_numbers, computer_capacity = computer_capacity)
    #corr.atoms["Vz"] = corr.atoms["Vz"] - 0.55
    #corr.show_density_plots()
    corr.compute_correlations()
    corr.result.to_pickle("~/res2.pkl")

    