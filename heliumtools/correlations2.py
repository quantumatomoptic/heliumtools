#!/usr/bin/env python
# -*- mode:Python; coding: utf-8 -*-

# ----------------------------------
# Created on the 5-4-2023 by Victor
#
# Developped by Victor, ...
#
# Last (big) change on the ... by ...
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

[1] https://pypi.org/project/fast-histogram/
[2] https://github.com/vaexio/vaex
[3] https://towardsdatascience.com/beyond-pandas-spark-dask-vaex-and-other-big-data-technologies-battling-head-to-head-a453a1f8cc13
"""


import numpy as np
import pandas as pd
import snoop

from tqdm import tqdm, trange
import copy
import time


class CorrelationHe2Style:
    """
    Classe CorrelationHe2Style. Prend en entrée un dataframe avec X, Y et T et calcule de corrélations etc...
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
        self.round_decimal = 7
        self.id = int(time.time())
        self.axis = ["Vx", "Vy", "Vz"]
        self.ref_frame_speed = {"Vx": 0, "Vy": 0, "Vz": 0}
        self.voxel_size = {"Vx": 0.1, "Vy": 0.1, "Vz": 0.1}
        self.voxel_numbers = {"Vx": 11, "Vy": 11, "Vz": 11}
        self.beams = {
            "A": {
                "Vx": {"size": 10, "position": 0},
                "Vy": {"size": 3, "position": 0},
                "Vz": {"size": 5, "position": 25},
            },
            "B": {
                "Vx": {"size": 10, "position": 0},
                "Vy": {"size": 3, "position": 0},
                "Vz": {"size": 5, "position": -25},
            },
        }
        self.__dict__.update(kwargs)
        # atoms : (mm,mm, ms) --> (mm/s, mm/s, mm/s)
        self.build_the_atoms_dataframe()
        self.cycles_array = atoms["Cycle"].unique()
        self.n_cycles = len(atoms["Cycle"].unique())
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

        self.atoms["X"] = 1000 * self.atoms["X"] / self.atoms["T"] + self.raman_kick
        self.atoms = self.atoms.rename(columns={"X": "Vx"})
        self.atoms["Y"] = 1000 * self.atoms["Y"] / self.atoms["T"]
        self.atoms = self.atoms.rename(columns={"Y": "Vy"})

        if (
            type(self.bec_arrival_time) == int
            or type(self.bec_arrival_time) == float
            or isinstance(self.bec_arrival_time, np.float64)
            or isinstance(self.bec_arrival_time, np.float32)
            or isinstance(self.bec_arrival_time, np.int32)
            or isinstance(self.bec_arrival_time, np.int64)
        ):
            l_fall = 0.5 * self.gravity * self.bec_arrival_time**2
            self.atoms["T"] = (
                0.5 * self.gravity * self.atoms["T"] - l_fall / self.atoms["T"]
            )
        elif isinstance(self.bec_arrival_time, pd.DataFrame):
            print("I change the bec arrival time for each cycle.")
            self.atoms = self.merge_dataframe_on_cycles(
                self.atoms, self.bec_arrival_time
            )
            l_fall = 0.5 * self.gravity * self.theoretical_arrival_time**2
            self.atoms["T"] = (
                0.5 * self.gravity * self.atoms["T"] - l_fall / self.atoms["T"]
            )
            # compute the speed of BECs
            self.atoms["BEC Arrival Time"] = (
                0.5 * self.gravity * self.atoms["BEC Arrival Time"]
                - l_fall / self.atoms["BEC Arrival Time"]
            )
            # take off the speed of each BEC
            self.atoms["T"] = (
                self.atoms["T"]
                - self.atoms["BEC Arrival Time"]
                + self.additionnal_bec_speed_ms
            )
            # drop the column with no more interest.
            self.atoms.drop("BEC Arrival Time", inplace=True, axis=1)

            ## On supprime ensuite les bad shots : là où il y a moins de 100 atomes
            df = self.bec_arrival_time[
                self.bec_arrival_time["Number of Atoms"] < self.bad_shot_limit
            ]
            for cycle, nb_atoms in zip(df["Cycle"], df["Number of Atoms"]):
                print(f"Delete cycle {cycle} with only {nb_atoms} atoms.")
                self.n_cycles -= 1

        else:
            print("###### /!\ Please modify build_the_atom_dataframe in correlation.py")
        self.atoms = self.atoms.rename(columns={"T": "Vz"})
        for axis in self.ref_frame_speed.columns:
            self.atoms[axis] -= self.ref_frame_speed[axis]

    def get_atoms_in_beams(self, beam):
        """
        Récupère  les atomes à l'intérieur de beam. Cela permet d'alléger les données à traiter.
        Si le beam est vide, la méthode ne fait rien. Le format du beam doit être {"Vx": {"position":120, "size":-120}}
        """
        if beam:
            for key, entry in beam.items():
                atoms_in_beam = self.atoms[
                    (
                        (self.atoms[key] <= entry["position"] + entry["size"] / 2.0)
                        & (self.atoms[key] > entry["position"] - entry["size"] / 2.0)
                    )
                ]
        return atoms_in_beam

    def update_atoms_in_beams(self):
        self.atomsA = self.get_atoms_in_beams(beam=self.beams["A"]).reset_index()
        self.atomsA["index"] = np.arange(0, len(self.atomsA))
        self.atomsB = self.get_atoms_in_beams(beam=self.beams["B"]).reset_index()
        self.atomsB["index"] = np.arange(0, len(self.atomsB))

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
        self.voxel_map_size = tuple([self.voxel_numbers[axis] for axis in self.axis])
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
        for column in ["G2AA", "G2BB", "G2AB"]:
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
        # atXY = atX.merge(atX, how="outer", on="Cycle")  # , on="index_A"), atABprime --> merge on cycles might be an issue when computing local correlation with uncorrelated atoms (normalization).
        atXY = atX.merge(atY, how="cross")
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
        atXY.drop(cols_to_remove, axis=1, inplace=True)
        G2XY, edges = np.histogramdd(
            atXY.to_numpy(), bins=self.voxel_map_size, range=self.voxel_map_range
        )
        return G2XY

    def compute_correlations(self):
        """Principal function of the code. L'idée du code est le suivant :
        * on parcourt les cycles de la séquence et pour chaque cycle, on récupère les atomes dans la ROI1 et la ROI2.
        * pour chaque cycle, on calcule la différence (corrélations locales) ou la somme (corrélation croisées) de l'ensemble des couple d'atomes.
        * on fait ensuite l'histogramme 3D de cet ensemble. On l'ajoute à l'histogramme total."""
        self.initialize_voxel_map_properties()
        self.update_atoms_in_beams()
        for cycle in tqdm(self.cycles_array):
            atA = self.atomsA[self.atomsA["Cycle"] == cycle]
            atB = self.atomsB[self.atomsB["Cycle"] == cycle]
            G2AA = self.get_G2(atA, atA, local=True)
            self.result["G2AA"] += G2AA.flatten()
            G2BB = self.get_G2(atB, atB, local=True)
            self.result["G2BB"] += G2BB.flatten()
            G2AB = self.get_G2(atA, atB, local=False)
            self.result["G2AB"] += G2AB.flatten()


if __name__ == "__main__":
    nat = 4
    ncycle = 3
    atoms = []
    cycle = []
    for ncy in range(ncycle):
        beamA = list(
            np.random.normal(loc=[0, 0, -25], scale=[10, 10, 1], size=(nat, 3))
        )
        beamB = list(np.random.normal(loc=[0, 0, 25], scale=[10, 10, 1], size=(nat, 3)))
        atoms += beamA + beamB
        cycle += list(ncy * np.ones(len(beamA) + len(beamB)))
    atoms = pd.DataFrame(atoms, columns=["Vx", "Vy", "Vz"])
    atoms["Cycle"] = cycle
    corr = CorrelationHe2Style(atoms)
    corr.compute_correlations()
