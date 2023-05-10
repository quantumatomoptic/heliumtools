#!/usr/bin/env python3
# -*- mode: Python; coding: utf-8 -*-

"""
@Author: victor
@Date:   20 April 2023 @ 19:01
@Last modified by:   victor
@Last modified time: 09 May 2023 @ 16:43

Comment :
"""


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
    """return the cycle from a path. If no cycle is found, return -1. Since April 2023, we save at each run of the experiment a unique id that give time in seconds since Helium1 Epoch time (aka first MOT with qcontrol3). 

    Parameters
    ----------
    path : str or path-like object vers un fichier .atoms
        chemin du .atoms
    """
    path = str(path)  # on s'assure que path soit bien un string (pas un Path object)
    try:
        json_path = path.replace(".atoms", ".json")
        with open(json_path) as f:
            list_of_metadatas = json.load(f)
        for element in list_of_metadatas:
            if element['name'] in ['seq', 'sequence', 'Sequence']:
                seq = element['value']
            elif element['name'] == 'cycle_id':
                cycle = element['value']
        return seq, cycle
    except:
        pass
    pattern = "[0-9][0-9][0-9]_[0-9][0-9][0-9][.]"
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
    selected_files = selected_files[0:N_files]
    print("Starting to gather atoms")
    for i in tqdm(range(N_files)):
        path = selected_files[i]
        X, Y, T = load_XYTTraw(path)
        seq, cycle = return_cycle_from_path(path)
        cycle_list = np.ones(len(X)) * cycle
        Xc = np.concatenate([Xc, X])
        Yc = np.concatenate([Yc, Y])
        Tc = np.concatenate([Tc, T])
        Cyclec = np.concatenate([Cyclec, cycle_list])
    # Je crée maintenant un dataframe
    data = np.transpose(np.array([Cyclec, Xc, Yc, Tc]))
    df = pd.DataFrame(data, columns=["Cycle", "X", "Y", "T"])
    return selected_files, df


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
    atom_files,
    ROI_for_fit={"T": {"min": 305.5, "max": 309.7}},
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
    def gaussian_function(x, mean, amplitude, standard_deviation):
        return amplitude * np.exp(-((x - mean) ** 2) / (2 * standard_deviation ** 2))
    list_of_arrival_time_max = [0 for i in atom_files]
    list_of_arrival_time= [0 for i in atom_files]
    number_of_atoms = [0 for i in atom_files]
    list_of_cycles = [0 for i in atom_files]
    print("Starting to gather arrival time of BECs")
    for i, path in enumerate(tqdm(atom_files)):
        X, Y, T = load_XYTTraw(path)
        number_of_atoms[i] = len(T)
        seq, cycle = return_cycle_from_path(path)
        list_of_cycles[i] = cycle
        if len(T) < 100:
            print(f"WARNING : shot {cycle} seems empty. I take it off the sequence.")
            list_of_arrival_time[i] = np.nan # 24/06/2022 & 17/05/2022
            list_of_arrival_time_max[i] = np.nan

        else:
            bin_heights, bin_borders = np.histogram(
                T,
                bins=np.arange(
                    ROI_for_fit["T"]["min"], ROI_for_fit["T"]["max"], histogramm_width
                ),
            )
            bin_centers = np.array(bin_borders[:-1] + np.diff(bin_borders) / 2)
            # find the position of the max

            list_of_arrival_time_max[i] = bin_centers[np.argmax(bin_heights)]
            p0 = [bin_centers[np.argmax(bin_heights)], np.max(bin_heights), np.std(T)]
            try:
                popt, pcov = curve_fit(gaussian_function, bin_centers, bin_heights, p0=p0)
                #perr = np.sqrt(np.diag(pcov))
            except:
                print(f"[WARN] : BEC not fitted (cycle {cycle} ; folder : {path}). Take the max as arrival time.")
                popt = p0
                #perr = [np.nan, np.nan, np.nan]
            list_of_arrival_time[i] = popt[0]

    df_arrival_time = pd.DataFrame(
        {
            "Cycle": list_of_cycles,
            "BEC Arrival Time": list_of_arrival_time,
            "Number of Atoms": number_of_atoms,
            "BEC Arrival Time with max":list_of_arrival_time_max,
            "BEC Arrival Time with fit":list_of_arrival_time
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
        ROI_for_fit=ROI_for_fit,
        histogramm_width=histogramm_width,
    )
    filename = os.path.join(directory, "arrival_times.pkl")
    df_arrival_time.to_pickle(filename)


def export_data_set_to_pickle(
    folder,
    ROI,
    find_arrival_times=False,
    n_max_cycles=1e8,
    histogramm_width=0.01,
    ROI_for_fit={"T": {"min": 306.2, "max": 309.7}},
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
    atom_files, atoms = load_atoms(folder, n_max_cycles=n_max_cycles)
    atoms_in_ROI = apply_ROI(atoms, ROI)
    filename_dataset = os.path.join(folder, "dataset.pkl")
    atoms_in_ROI.to_pickle(filename_dataset)
    df_parameters = gather_saved_sequence_parameters(folder)
    if find_arrival_times:
        df_arrival_times = obtain_arrival_times(
            atom_files, histogramm_width=histogramm_width, ROI_for_fit=ROI_for_fit
        )
        df_arrival_times = pd.merge(df_arrival_times, df_parameters, on="Cycle")
        filename = os.path.join(folder, "arrival_times.pkl")
        df_arrival_times.to_pickle(filename)

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
            if element["name"] not in ["sequence", "cycle", "cycle_id", "Sequence", "Cycle"]:
                column_names.append("{} ({})".format(element["name"], element["unit"]))
                column_values.append(element["value"])
        new_df = pd.DataFrame([column_values], columns=column_names)
        dataframe = pd.concat([dataframe, new_df])
    dataframe.reset_index(drop=True, inplace=True)
    filename_parameters = os.path.join(folder, "parameters.pkl")
    dataframe.to_pickle(filename_parameters)
    return dataframe


if __name__ == "__main__":
    folder = Path("/mnt/manip_E/2022/11/18/083")
    export_data_set_to_pickle(
        folder,
        ROI={"T": {"min": 306.2, "max": 309.7}},
        find_arrival_times=True,
        n_max_cycles=3,
    )
