# !/usr/bin/env python
#  -*-mode:Python; coding:utf8 -*-

# ------------------------------------
# Created on Mo Oct 2023 by Victor Gondret
# Contact : victor.gondret@institutoptique
#
# MIT Copyright (c) 2023 - Helium1@LCF
# Institut d'Optique Graduate School
# Université Paris-Saclay
# ------------------------------------
#
"""
Decription of bell_database.py 

Please document your code ;-).
"""


import os, glob, re
import pandas as pd
from heliumtools.misc.gather_data import apply_ROI
import logging
import math

log = logging.getLogger(__name__)  #


def data_filter(data, bec_arrival_times, filters):
    selec_bec_arrival_times = apply_ROI(bec_arrival_times, filters)
    selected_data = data[data["Cycle"].isin(selec_bec_arrival_times["Cycle"])]
    return selected_data, selec_bec_arrival_times


def add_sequence_to_bell_folder(bell_folder, sequence_folder, do_export=False):
    """
    Function that add a sequence to the bell folder database.  It contains pkl file with names
    atoms_theta_{angle}_degrees.pkl containing all atoms with the corresponding angle
    and the corresponding file bec_arr_time_theta_{angle}_degrees.pkl with BEC arrival time.
    Parameter
    ---------
    sequence : string or path to load the sequence.
    """
    if do_export:
        export_data_set_to_pickle(
            sequence_folder,
            find_arrival_times=True,
            histogramm_width=0.01,
            width_saturation=0.1,
            ROI={"T": {"min": 312, "max": 322}},
            ROI_for_fit={"T": {"min": 305, "max": 310}},
            supplementary_rois=[{"T": [450, 475]}],
            metadata=["picoscope"],
        )
    try:
        atoms = pd.read_pickle(os.path.join(sequence_folder, "dataset.pkl"))
        bec_arrival_times = pd.read_pickle(
            os.path.join(sequence_folder, "arrival_times.pkl")
        )
    except Exception as err:
        msg = f"Atoms and arrival time loading failed from folder {sequence_folder}. "
        msg += f"Please make sure they do exist. Error descirption : {err}"
        log.error(msg)
        return
    sequence_id = sequence_folder
    bec_arrival_times["Sequence Id"] = sequence_id
    sequence_angles = bec_arrival_times["thetAspect (dimensionless)"].unique()
    for angle in sequence_angles:
        if math.isnan(angle):
            continue
        filename_suffix = "{:0>4d}".format(round(angle))
        atoms_filename = os.path.join(
            bell_folder, "atoms_theta_{}_degrees.pkl".format(filename_suffix)
        )
        bec_filename = os.path.join(
            bell_folder, "bec_arr_time_theta_{}_degrees.pkl".format(filename_suffix)
        )
        if os.path.exists(atoms_filename):
            atoms_angle_fixed = pd.read_pickle(atoms_filename)
        else:
            atoms_angle_fixed = pd.DataFrame({"Cycle": [], "X": [], "Y": [], "T": []})
        if os.path.exists(bec_filename):
            bec_angle_fixed = pd.read_pickle(bec_filename)
        else:
            bec_angle_fixed = pd.DataFrame({"Sequence Id": [], "Cycle": []})
        # We want to add data to
        if sequence_id in bec_angle_fixed["Sequence Id"].unique():
            print(
                f"[WARNING] The sequence {sequence_id} you want to add was already in "
                + f"the bell Database at {round(angle)} us. I suppress all atoms corresponding"
                + " to this sequence and add those new ones."
            )
            bec_angle_fixed = bec_angle_fixed[
                ~(bec_angle_fixed["Sequence Id"] == sequence_id)
            ]
            atoms_angle_fixed = atoms_angle_fixed[
                atoms_angle_fixed["Cycle"].isin(bec_angle_fixed["Cycle"])
            ]
        bec_to_add = bec_arrival_times[
            bec_arrival_times["thetAspect (dimensionless)"] == angle
        ]
        atoms_to_add = atoms[atoms["Cycle"].isin(bec_to_add["Cycle"])]
        bec_angle_fixed = pd.concat([bec_angle_fixed, bec_to_add])
        atoms_angle_fixed = pd.concat([atoms_angle_fixed, atoms_to_add])
        bec_angle_fixed.to_pickle(bec_filename)
        atoms_angle_fixed.to_pickle(atoms_filename)


def load_all_bell_database(bell_folder):
    bec_files = glob.glob(bell_folder + "/bec_arr_time_theta_*.pkl")
    atoms_files = glob.glob(bell_folder + "/atoms_theta_*.pkl")
    raw_atoms = []
    raw_arrival_times = []
    for bec_file, atom_file in zip(bec_files, atoms_files):
        raw_atoms.append(pd.read_pickle(atom_file))
        raw_arrival_times.append(pd.read_pickle(bec_file))
    atoms = pd.concat(raw_atoms).reset_index(drop=True)
    bec_arrival_times = pd.concat(raw_arrival_times).reset_index(drop=True)
    return atoms, bec_arrival_times


def load_all_bec_arrival_times(bell_folder):
    bec_files = glob.glob(bell_folder + "/bec_arr_time_theta_*.pkl")
    raw_arrival_times = []
    for bec_file in bec_files:
        raw_arrival_times.append(pd.read_pickle(bec_file))
    return pd.concat(raw_arrival_times).reset_index(drop=True)


def get_all_angle(bell_folder):
    bec_files = glob.glob(bell_folder + "/bec_arr_time_theta_*.pkl")
    angle_list = []
    for filename in bec_files:
        pattern = "[0-9][0-9][0-9][0-9]_degrees.pkl"
        angle_list.append(
            int(re.findall(pattern, filename)[0].replace("_degrees.pkl", ""))
        )
    return angle_list


def load_data_angle(bell_folder, angle):
    bec_file = os.path.join(
        bell_folder + "/bec_arr_time_theta_{:0>4d}_degrees.pkl".format(angle)
    )
    atoms_file = os.path.join(
        bell_folder + "/atoms_theta_{:0>4d}_degrees.pkl".format(angle)
    )
    if os.path.exists(bec_file) and os.path.exists(atoms_file):
        atoms = pd.read_pickle(atoms_file).reset_index(drop=True)
        bec_arrival_times = pd.read_pickle(bec_file).reset_index(drop=True)
        return atoms, bec_arrival_times
    else:
        print("[WARNING]: No data were found with angle of {} us.".format(angle))
        return pd.DataFrame(), pd.DataFrame()
