#!/usr/bin/env python
# -*- mode:Python; coding: utf-8 -*-

# ----------------------------------
# Created on the 20-6-2023 by Victor, Rui
#
# Copyright (c) 2022 - Helium1@LCF
# ----------------------------------


"""
Content of data_builder.py
-----------------------------

"""

import numbers
import numpy as np
import pandas as pd
from tqdm import tqdm, trange
import copy, random
import logging

from heliumtools.misc.gather_data import apply_ROI, apply_ROD


class DataBuilder:
    """
    Classe DataBuilder. Takes as an input a pandas dataframe with 4 columns [Cycle, X, Y et T] with X,Y in mm and T in ms
    and builds a dataframe in mm/s for all speeds. The name os this dataframe is self.atoms.
    In the computation of the atoms velocity it is assumed that the Raman is applied right after (a few us) the switch off of the dipole trap.

    Mandatory parameters
    --------------------
    atoms : pandas dataframe
        Must have 4 columns in which each columns name is "Cycle" ; "X" ; "Y" et "T".

    Other attributs
    --------------------

    ROI : dict
        Define a Region Of Interest to select a part of the DataFrame. The entries of this dictionary must be in [Cycle, Vx, Vy, Vz, Vperp, theta, Vrho, phi]
        The apply_ROI method is applied immediately after the dataframe atom is built.

    ROD : dictionnaire, Region Of Desinterest
        To select only atoms outside the ROD. Same remark as for ROI and same format.

    bec_arrival_time : float or pandas dataframe
        if a float, it should be the condensate arrival time in ms. Else, it should be a dataframe with at least two collumns ["Cycle","BEC Arrival Time"]. 
        In total it can have the collumns ["Cycle","BEC Center X","BEC Center X","BEC Arrival Time"] where BEC Arrival Time is the arrival time in ms at each cycle 
        and BEC Center X and BEC Center X are the positions in mm at each cycle.

    raman_kick : dict
        Raman kick in mm/s in each direction

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
        # condensate arrival time in ms
        self.bec_arrival_time = 307.763  
        # Raman kick in mm/s. It can also be the initial velocity
        self.raman_kick = {"Vx": -42.5 , "Vy" : 0.0 , "Vz" : -4.55}  
        # gravitational acceleration in m/s^2 
        self.gravity = 9.81 
        # velocity of reference frame in mm/s
        self.refFrame = {"Vx": 0, "Vy": 0, "Vz": 0}
        # define region of interest to keep
        self.ROI = {}
        # define region of interest to remove
        self.ROD = {}
        # True or false you you want to recenter dataframe with respect to BEC
        self.recenter_with_bec_arrival_time = {"Vx": True, "Vy": True, "Vz": True}
        # if there is a copy of atoms dataframe. For bootstrap
        self.is_there_a_copy_of_atoms = False

        # update arguments of class
        self.__dict__.update(kwargs)
        # build compute dataframe with atoms velocities
        self.build_the_atoms_dataframe()

    def build_the_atoms_dataframe(self):
        """
        This method builds the dataframe containing all the positions of the atoms. It builds the dataframe self.atoms, whose speeds are expressed in mm/s from the initial dataframe.
        """
        # Clean the self.atoms dataframe : should contain only 4 columns (to avoid mixing when merging).
        for column in self.atoms.columns:
            if column not in ["X", "Y", "T", "Cycle"]:
                self.atoms.drop(column, inplace=True, axis=1)

        # update units of gravity to mm/ms^2
        gravity = self.gravity*1e-3
        # update units of reference frame and raman kick to mm/ms
        raman_kick = dict()
        refFrame = dict()
        for i in ["Vx","Vy","Vz"]:
            raman_kick[i] = self.raman_kick[i]*1e-3
            refFrame[i] = self.refFrame[i]*1e-3
        
        # if self.bec_arrival_time is a number
        if isinstance(self.bec_arrival_time, numbers.Number):
            # compute initial position of the atoms using BEC arrival time 
            l_fall = -raman_kick["Vz"]*self.bec_arrival_time + 0.5*gravity*self.bec_arrival_time**2

            # compute velocity along x-axis
            self.atoms["X"] = self.atoms["X"]/self.atoms["T"] - raman_kick["Vx"]
            self.atoms = self.atoms.rename(columns={"X": "Vx"})

            # compute velocity along y-axis
            self.atoms["Y"] = self.atoms["Y"]/self.atoms["T"] - raman_kick["Vy"]
            self.atoms = self.atoms.rename(columns={"Y": "Vy"})

            # compute velocity along z-axis
            self.atoms["T"] =  0.5*gravity*self.atoms["T"] - l_fall/self.atoms["T"] - raman_kick["Vz"]
            self.atoms = self.atoms.rename(columns={"T": "Vz"})

            # update units to mm/s and apply reference frame
            for axis in ["Vx","Vy","Vz"]:
                self.atoms[axis] = self.atoms[axis]*1e3
                if axis in self.ref_frame_speed.keys():
                    self.atoms[axis] = self.atoms[axis] - self.ref_frame_speed[axis]
                else:
                    self.ref_frame_speed[axis] = 0.0

        elif isinstance(self.bec_arrival_time, pd.DataFrame):
            if "Cycle" in self.bec_arrival_time.keys() and "BEC Arrival Time" in self.bec_arrival_time.keys():

                # compute initial position of the atoms using BEC arrival time 
                self.bec_arrival_time["l_fall"] = -raman_kick["Vz"]*self.bec_arrival_time["BEC Arrival Time"]+ 0.5*gravity*np.power(self.bec_arrival_time["BEC Arrival Time"],2)

                # merge self.bec_arrival_time in self.atoms
                self.atoms = self.atoms.merge(self.bec_arrival_time,on = "Cycle")

                # compute velocity along x-axis
                self.atoms["X"] = self.atoms["X"]/self.atoms["T"] - raman_kick["Vx"]
                self.atoms = self.atoms.rename(columns={"X": "Vx"})

                # compute velocity along y-axis
                self.atoms["Y"] = self.atoms["Y"]/self.atoms["T"] - raman_kick["Vy"]
                self.atoms = self.atoms.rename(columns={"Y": "Vy"})

                # compute velocity along z-axis
                self.atoms["T"] = 0.5*gravity*self.atoms["T"] - self.atoms["l_fall"]/self.atoms["T"] - raman_kick["Vz"]
                self.atoms = self.atoms.rename(columns={"T": "Vz"})

                # compute BEC velocity along x-axis
                if "BEC Center X" in self.bec_arrival_time.columns:
                    self.atoms["BEC Center X"] = self.atoms["BEC Center X"]/self.atoms["BEC Arrival Time"] - raman_kick["Vx"]
                else:
                    self.atoms["BEC Center X"] = 0.0

                # compute BEC velocity along y-axis 
                if "BEC Center Y" in self.bec_arrival_time.columns:
                    self.atoms["BEC Center Y"] = self.atoms["BEC Center Y"]/self.atoms["BEC Arrival Time"] - raman_kick["Vy"]                    
                else:
                    self.atoms["BEC Center Y"] = 0.0

                #  compute BEC velocity along z-axis 
                self.atoms["BEC Center Z"] = 0.5*gravity*self.atoms["BEC Arrival Time"] - self.atoms["l_fall"]/self.atoms["BEC Arrival Time"] - raman_kick["Vz"]
        
                # Recenter velocities using BEC velocity
                for Vj, AX in zip(["Vx", "Vy", "Vz"], ["X", "Y", "Z"]):
                    if self.recenter_with_bec_arrival_time[Vj]:
                        self.atoms[Vj] = self.atoms[Vj] - self.atoms["BEC Center " + AX]

                # drop columns with no more interest (for later calculations)
                for column in self.atoms.columns:
                    if column not in ["Vz", "Vx", "Vy", "Cycle", "Vperp", "theta"]:
                        self.atoms.drop(column, inplace=True, axis=1)

                # update units to mm/s and apply reference frame
                for axis in ["Vx","Vy","Vz"]:
                    self.atoms[axis] = self.atoms[axis]*1e3
                    if axis in self.ref_frame_speed.keys():
                        self.atoms[axis] = self.atoms[axis] - self.ref_frame_speed[axis]
                    else:
                        self.ref_frame_speed[axis] = 0.0

            else:
                logging.error("[ERROR] From build_the_atoms_dataframe : the bec_arrival_time dataframe is not valid.")
    
        else:
            logging.error("[ERROR] From build_the_atoms_dataframe : the bec_arrival_time instance is not recognized.")


        self.cycles_array = self.atoms["Cycle"].drop_duplicates()
        self.n_cycles = len(self.cycles_array)

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
            print("[Warning] : I just saved a copy of the atom dataframe because you will destruct your original dataframe.")
        new_atoms = []
        for n in range(self.n_cycles):
            cycle = random.choice(self.cycles_array_copy)
            df = copy.deepcopy(self.atoms_dataframe_copy[self.atoms_dataframe_copy["Cycle"] == cycle])
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

    def rotate_inertial_frame(self, anglexy=0, anglexz=0, angleyz=0):
        """
        Rotate the atoms dataframe in 3D space using rotation matrices.

        Parameters
        ----------
        anglexy : float
            Clockwise rotation angle in the XY plane, in degrees.
        anglexz : float
            Clockwise rotation angle in the XZ plane, in degrees.
        angleyz : float
            Clockwise rotation angle in the YZ plane, in degrees.
        """
        # Convert degrees to radians for rotation calculations
        anglexy_rad = np.radians(anglexy)
        anglexz_rad = np.radians(anglexz)
        angleyz_rad = np.radians(angleyz)

        # Extract original velocity components
        X = self.atoms["Vx"]
        Y = self.atoms["Vy"]
        Z = self.atoms["Vz"]  # Assuming you also have a Z component

        # Rotation in XY Plane
        cos_xy = np.cos(anglexy_rad)
        sin_xy = np.sin(anglexy_rad)
        X_new = X * cos_xy + Y * sin_xy
        Y_new = -X * sin_xy + Y * cos_xy

        # Rotation in XZ Plane
        cos_xz = np.cos(anglexz_rad)
        sin_xz = np.sin(anglexz_rad)
        Z_new = Z * cos_xz - X_new * sin_xz
        X_new = X_new * cos_xz + Z * sin_xz

        # Rotation in YZ Plane
        cos_yz = np.cos(angleyz_rad)
        sin_yz = np.sin(angleyz_rad)
        Y_new = Y_new * cos_yz + Z_new * sin_yz
        Z_new = -Y_new * sin_yz + Z_new * cos_yz

        # Update the atoms DataFrame with new velocity components
        self.atoms["Vx"] = X_new
        self.atoms["Vy"] = Y_new
        self.atoms["Vz"] = Z_new
