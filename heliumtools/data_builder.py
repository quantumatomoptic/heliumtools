#!/usr/bin/env python
# -*- mode:Python; coding: utf-8 -*-

# ----------------------------------
# Created on the 20-6-2023 by Victor
#
# Copyright (c) 2022 - Helium1@LCF
# ----------------------------------
#
"""
Content of data_builder.py
-----------------------------

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm, trange
import copy, random
import matplotlib.cm as cm
import matplotlib.colors as colors
import time
import logging
from .misc.gather_data import apply_ROI, apply_ROD


class DataBuilder:
    """
    Classe DataBuilder. Takes as an input a pandas dataframe with 4 columns Cycle, X, Y et T and built a dataframe in mm/s for all speeds.

    Mandatory parameters
    --------------------
    atoms : pandas dataframe
        Must have 4 columns in which each columns name is "Cycle" ; "X" ; "Y" et "T".



    Other attributs
    --------------------

    ROI : dict
        Define a Region Of Interest to select a part of the DataFrame. Note that any entry of this dictionary must be in the builted atom dataframe (must be among Cycle, Vx, Vy, Vz, Vperp, theta, Vrho, phi)
        Note that the apply_ROI method is applied immediately after the dataframe atom is built.

    ROD : dictionnaire, Region Of Desinterest
        pour ne sélectionner que les atomes en dehors de la ROD. Même remarque que pour ROI et même format.

    bec_arrival_time : float, temps d'arrivée du condensat en ms

    raman_kick : float, kick raman en mm/s

    Some Methods
    --------------------
    apply_ROI() / apply_ROD() : select atoms only in ROI / outside ROD.


    show_density : plot the density of the atoms dataframe

    compute_spherical_coordinates

    """

    def __init__(self, atoms, **kwargs):
        """
        Object initialization, sets parameters as the user defined, build the atoms dataframe and apply ROD and ROI.
        """
        self.atoms = copy.deepcopy(atoms)
        self.cycles_array = self.atoms["Cycle"].unique()
        self.n_cycles = len(self.atoms["Cycle"].unique())
        self.bec_arrival_time = 307.763  # temps d'arrivée du BEC, en ms
        self.theoretical_arrival_time = 307.763  # 24/06/2022 & 17/05/2022
        self.raman_kick = 42.5  # mm/s, kick Raman
        self.gravity = 9.81
        self.bad_shot_limit = 100
        self.ref_frame_speed = {"Vx": 0, "Vy": 0, "Vz": 0}
        self.ROI = {}  # Region Of Interest
        self.ROD = {}  # Region Of Desinterest.
        self.recenter_with_bec_arrival_time = {"Vx": True, "Vy": True, "Vz": True}
        self.construct_atoms = True
        self.is_there_a_copy_of_atoms = False
        self.bec_arrival_time = copy.deepcopy(self.bec_arrival_time)
        self.id = int(time.time())
        self.__dict__.update(kwargs)
        if self.construct_atoms:
            self.build_the_atoms_dataframe()
        self.apply_ROI()  # Keep only atoms in ROI
        self.apply_ROD()  # Take off atoms in Region Of Desinterest.
        self.cycles_array = self.atoms["Cycle"].unique()
        self.n_cycles = len(self.atoms["Cycle"].unique())

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

            ## --- Merge dataframe
            # clean before merging
            # drop columns with no more interest (for later calculations)
            bec_arr_times = copy.deepcopy(self.bec_arrival_time)
            for column in bec_arr_times:
                if column not in [
                    "BEC X",
                    "BEC Y",
                    "BEC Z",
                    "Cycle",
                ]:
                    bec_arr_times.drop(column, inplace=True, axis=1)

            self.atoms = self.merge_dataframe_on_cycles(self.atoms, bec_arr_times)
            # take off the speed of each BEC : recenter
            for Vj, AX in zip(["Vx", "Vy", "Vz"], ["X", "Y", "Z"]):
                if self.recenter_with_bec_arrival_time[Vj]:
                    self.atoms[Vj] = self.atoms[Vj] - self.atoms["BEC " + AX]
                elif Vj == "Vx":  # raman kick to be taken into account
                    self.atoms[Vj] = self.atoms[Vj] + self.raman_kick
            # drop columns with no more interest (for later calculations)
            for column in self.atoms.columns:
                if column not in ["Vz", "Vx", "Vy", "Cycle", "Vperp", "theta"]:
                    self.atoms.drop(column, inplace=True, axis=1)

        else:

            logging.error(
                "[ERROR] From build_the_atoms_dataframe : the bec_arrival_time instance is not recognized."
            )

        for axis in ["Vx", "Vy", "Vz"]:
            if axis in self.ref_frame_speed:
                # logging.info(
                #    f"[INFO] : Reference frame is moving at {self.ref_frame_speed[axis]} mm/s along the {axis} axis."
                # )
                self.atoms[axis] -= self.ref_frame_speed[axis]
            else:
                self.ref_frame_speed[axis] = 0

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

    def compute_spherical_coordinates(self):
        """Compute spherical coordinates for the atom dataframe. Note that the call of that function might overlap with cylindrical coordinates."""
        self.atoms["rho"] = np.sqrt(
            self.atoms["Vx"] ** 2 + self.atoms["Vy"] ** 2 + self.atoms["Vz"] ** 2
        )
        self.atoms["theta"] = np.arccos(
            self.atoms["Vz"] / self.atoms["rho"]
        )  # [0, pi] range

        self.atoms["phi"] = np.arctan2(self.atoms["Vy"], self.atoms["Vx"]) % (2 * np.pi)
        # phi entre -pi et pi --> je mets entre 0 et 2pi

    def update_referential_speed(self, new_referential_speed: dict):
        """This methods update the speed of the inertial frame taking into account the old inertial frame. This means that this frame is absolute with respect to the detected speed of atoms.

        Parameters
        ----------
        new_referential_speed : dict
            dictionary with entries whos element are Vx, Vy and/or Vz and a float.
        """
        for axis in ["Vx", "Vy", "Vz"]:
            if axis not in new_referential_speed:
                new_referential_speed[axis] = 0
            self.atoms[axis] -= new_referential_speed[axis] - self.ref_frame_speed[axis]

        self.ref_frame_speed = copy.deepcopy(new_referential_speed)
        self.compute_cylindrical_coordinates()

    def return_dictionary_correlation_property(self) -> dict:
        """Return a dictionary with all the parameter of the simulation.

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

    def apply_ROI(self):
        """
        Modifie le dataframe "atoms" en appliquant la ROI. Cela permet d'alléger les données à traiter.
        Si la ROI est vide, la méthode ne fait rien. Le format de la ROI doit être {"Vx": {"max":120, "min":-120}}
        """
        self.atoms = apply_ROI(self.atoms, self.ROI)

    def apply_ROD(self):
        """
        Modifie le dataframe "atoms" en appliquant la region of desinterest i.e. en sélectionnant les atomes autre que ceux dans la ROD. Si la ROD est vide, la méthode ne fait rien.
        """
        self.atoms = apply_ROD(self.atoms, self.ROD)

    def rotate_inertial_frame(self, angle: float):
        """rotate the atoms dataframe in the (Vx, Vy) plane with an angle.

        Parameters
        ----------
        angle : float
            rotation angle, in degree
        """
        X = self.atoms["Vx"]
        Y = self.atoms["Vy"]
        self.atoms["Vx"] = X * np.cos(2 * np.pi * angle / 360) + Y * np.sin(
            2 * np.pi * angle / 360
        )
        self.atoms["Vy"] = -Y * np.sin(2 * np.pi * angle / 360) + Y * np.cos(
            2 * np.pi * angle / 360
        )
