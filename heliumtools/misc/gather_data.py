#!/usr/bin/env python
# -*- mode:Python; coding: utf-8 -*-

# ----------------------------------
# Created on the Tue Jun 28 2022
# Author : Victor
# Copyright (c) 2022 - Helium1@LCF
# ----------------------------------
#
"""
Content of gather_data.py

In this file are defined some functions to gather data.
"""
import glob, os, re, json
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
from tqdm import tqdm


def select_atoms_in_folder(folder):
    """Renvoie l'ensemble des fichiers .atoms de folder.

    Parameters
    ----------
    folder : path like object, NOT STRING
        chemin dans lequel aller chercher les .atoms

    Returns
    -------
    list of paths
        liste des chemins avec des .atoms dans le fichier.
    """
    folder = Path(folder)
    ext = [".atoms"]
    path_list = sorted(
        [
            path.as_posix()
            for path in filter(lambda path: path.suffix in ext, folder.glob("*"))
        ]
    )  # string object
    return path_list


def return_cycle_from_path(path):
    """return the cycle from a path. If no cycle is found, return -1

    Parameters
    ----------
    path : _type_
        _description_
    """
    pattern = "[0-9][0-9][0-9]_[0-9][0-9][0-9][.]"
    path = str(path)  # on s'assure que path soit bien un string (pas un Path object)
    result = re.findall(pattern, path)
    if len(result) == 0:
        pattern = "[0-9][0-9][0-9]_[0-9][0-9][0-9][0-9][.]"
        # On test si jamais le cycle est supérieur à 1000.
        result = re.findall(pattern, path)
        if len(result) == 0:
            return -1
    if len(result) > 1:
        print("Strange, I found two regular expression : {}".format(result))
        expr = result[-1]
    else:
        expr = result[0]
    # On a maitnenant expr qui est de la forme '003_002.' --> on récupère juste 002 dans ce cas
    expr = expr.replace(".", "")
    seq = int(expr[0:3])
    expr = expr[4:]
    cycle = int(expr)
    return seq, cycle


def load_atoms(folder, n_max_cycles=1e8):
    """Charge tous les .atoms.

    Parameters
    ----------
    folder : path like
        chemin vers le dossier contenant tous les .atoms
    """
    selected_files = select_atoms_in_folder(Path(folder))
    N_files = min([len(selected_files), n_max_cycles])
    Xc, Yc, Tc, Cyclec = [], [], [], []
    print("Starting to gather atoms")
    for i in tqdm(range(N_files)):
        path = selected_files[i]
        X, Y, T = load_XYTTraw(path)
        cycle = np.ones(len(X)) * (1 + i)
        Xc = np.concatenate([Xc, X])
        Yc = np.concatenate([Yc, Y])
        Tc = np.concatenate([Tc, T])
        Cyclec = np.concatenate([Cyclec, cycle])
    # Je crée maintenant un dataframe
    data = np.transpose(np.array([Cyclec, Xc, Yc, Tc]))
    df = pd.DataFrame(data, columns=["Cycle", "X", "Y", "T"])
    return N_files, df


def load_XYTTraw(path):
    """Retourne les positions X, Y, T d'un cycle d'une séquence.

    Parameters
    ----------
    path : str
        nom du fichier SANS EXTENSION

    Returns
    -------
    X, Y, T numpy arrays of same size
    """
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


def apply_ROI(atoms, ROI):
    """
    Modifie le dataframe "atoms" en appliquant la ROI. Cela permet d'alléger les données à traiter.
    Si la ROI est vide, la méthode ne fait rien.

    Parameters
    ----------
    atoms : pandas dataframe
        pandas dataframe avec les
    ROI : dic
        dictionnaire du type {"T": {"max":120, "min":-120}} où les entrées correspondent aux entrées du dataframe atoms
    """
    if ROI:
        for key, entry in ROI.items():
            atoms = atoms[((atoms[key] <= entry["max"]) & (atoms[key] > entry["min"]))]
    return atoms


def obtain_arrival_times(
    atoms,
    histogramm_width=0.01,
):
    """Pour chaque cycle du dataframe atoms, on fait un histogramme des temps d'arrivée des atomes puis on considère que le max de cet histogramme correspond au temps d'arrivé du BEC. Le paramètre histogramm_width permet d'ajuster la largeur des pics de l'histogramme (ms)
    WARNING : cette fonction est longue et non optimisée.

    Parameters
    ----------
    atoms : list of string path
        liste des cycles avec les chemins allant au .atoms.
    ROI : dic
        Region d'interet pour récupérer le temps d'arrivée du BEC
    histogramm_width : float (ms)
        largeur de chaque pic de l'histogramme.

    Return
    ------
    dataframe containing BEC arrival time & # of atoms
    """
    list_of_cycles = atoms["Cycle"].unique()
    list_of_arrival_time = [0 for i in list_of_cycles]
    number_of_atoms = [0 for i in list_of_cycles]
    print("Starting to gather arrival time of BECs")
    for i, cycle in enumerate(tqdm(list_of_cycles)):
        df = atoms[atoms["Cycle"] == cycle]
        if len(df) < 100:
            print(f"WARNING : shot {cycle} seems empty. I take it off the sequence.")
            list_of_arrival_time[i] = 308.07  # 24/06/2022 & 17/05/2022
            number_of_atoms[i] = len(df["T"])
        else:
            bin_heights, bin_borders, _ = plt.hist(
                df["T"],
                bins=np.arange(np.min(df["T"]), np.max(df["T"]), histogramm_width),
            )
            plt.close()
            bin_centers = np.array(bin_borders[:-1] + np.diff(bin_borders) / 2)
            # find the position of the max
            bec_arrival_time = bin_centers[np.argmax(bin_heights)]
            list_of_arrival_time[i] = bec_arrival_time
            number_of_atoms[i] = len(df["T"])
    df_arrival_time = pd.DataFrame(
        {
            "Cycle": list_of_cycles,
            "BEC Arrival Time": list_of_arrival_time,
            "Number of Atoms": number_of_atoms,
        }
    )
    return df_arrival_time


def function_find_arrival_times(
    atoms,
    directory,
    ROI_for_fit={"T": {"min": 306.2, "max": 309.7}},
    histogramm_width=0.01,
):
    """_summary_

    Parameters
    ----------
    atoms : pandas dataframe.
        dataframe of atoms, column T, X, Y et Cycle.
    directory : string or pathlike object
        directory in which one wants to save the pickle file
    ROI_for_fit : dict, optional
        ROI pour fitter les temps d'arrivée, by default {"T": {"min": 307, "max": 309}}
    histogramm_width : float (ms), optional
        largeur des pics de l'histogramme, en ms, by default 0.01
    """
    # STEP 2 : find arrival times and save them
    atoms_for_arrival_time = apply_ROI(atoms, ROI_for_fit)
    df_arrival_time = obtain_arrival_times(
        atoms_for_arrival_time,
        histogramm_width=histogramm_width,
    )
    filename = os.path.join(directory, "arrival_times.pkl")
    df_arrival_time.to_pickle(filename)


def export_data_set_to_pickle(
    folder, ROI, find_arrival_times=False, n_max_cycles=1e8, histogramm_width=0.01
):
    """Exporte le dataset folder comme pickle.

    Parameters
    ----------
    folder : path like
        chemin vers le dossier contenant tous les .atoms
    ROI : dictionnaire
     Exemple : {"T": {"min": 300, "max": 350}}

    """
    ### STEP 1 : gather data and save it
    N_files, atoms = load_atoms(folder, n_max_cycles=n_max_cycles)
    atoms_in_ROI = apply_ROI(atoms, ROI)
    filename_dataset = os.path.join(folder, "dataset.pkl")
    atoms_in_ROI.to_pickle(filename_dataset)

    if find_arrival_times:
        df_arrival_times = obtain_arrival_times(atoms, histogramm_width)
        df_parameters = gather_saved_sequence_parameters(folder)
    # df_arrival_time = pd.merge(df_arrival_times, df_parameters, on="Cycle")
    filename = os.path.join(directory, "arrival_times.pkl")
    df_arrival_time.to_pickle(filename)
    return filename_dataset


def gather_saved_sequence_parameters(folder):
    """Cette fonction récupère tous les paramètres sauvegardés pour chaque cycle de
    la séquence et les mets dans un dataframe avec le numéro du cycle.

    Parameters
    ----------
    folder : path like
        chemin vers le dossier contenant tous les .atoms
    """
    atoms = select_atoms_in_folder(folder)
    dataframe = pd.DataFrame()
    for atom_name in atoms:
        filename = atom_name.replace(".atoms", ".json")
        f = open(filename)
        data = json.load(f)
        column_names = ["Sequence", "Cycle"]
        seq, cycle = return_cycle_from_path(atom_name)
        column_values = [seq, cycle]
        for element in data:
            column_names.append("{} ({})".format(element["name"], element["unit"]))
            column_values.append(element["value"])
        new_df = pd.DataFrame([column_values], columns=column_names)
        dataframe = pd.concat([dataframe, new_df])
    dataframe.reset_index(drop=True, inplace=True)
    filename_parameters = os.path.join(folder, "parameters.pkl")
    dataframe.to_pickle(filename_parameters)
    return dataframe


if __name__ == "__main__":
    folder = Path("/mnt/manip_E/2022/05/18/089")
    export_data_set_to_pickle(folder, find_arrival_times=True)
