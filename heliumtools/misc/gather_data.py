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
# Authors : Charlie, Victor
# Copyright (c) 2022 - Helium1@LCF
# ----------------------------------
#
"""
Content of gather_data.py

In this file are defined some functions to gather data.
"""
from scipy.optimize import curve_fit
import glob, os, re, json, random
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
from tqdm import tqdm


def gaussian_function(x, mean, amplitude, standard_deviation, offset):
    return offset + amplitude * np.exp(
        -((x - mean) ** 2) / (2 * standard_deviation**2)
    )


def apply_roi(df, roi):
    for key, value in roi.items():
        df = df[df[key] > value["min"]]
        df = df[df[key] < value["max"]]
    return df


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
        if ".atoms" in path:
            json_path = path.replace(".atoms", ".json")
        elif ".json" in path:
            json_path = path
        with open(json_path) as f:
            list_of_metadatas = json.load(f)
        for element in list_of_metadatas:
            if element["name"] in ["seq", "sequence", "Sequence"]:
                seq = element["value"]
            elif element["name"] == "cycle_id":
                cycle = element["value"]
                run_id = element["value"]
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
    return (
        seq,
        cycle,
    )


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


def apply_ROD(df, ROD):
    """This function returns athe atoms dataframe such that all elements are OUTSIDE the ROD dictionary range. If the ROD is an empty dictionary, this function returns the initial dataframe.

    Parameters
    ----------
    df : pandas dataframe
        DataFrame with all atoms.
    ROI : dic
        dictionary for which every entry (for exemple 'T') matches a column of the df dataframe. The function returns a dictionary with the same number of columns as df but for which every line is NOT in a range required by the ROD dictionary, i.e. a maximum and a minimum value. Since a recent update, each entry of the dictionary can be a tuple or a list with two number or a dictionary with entries "min/max" or 'range'.

    Returns
    -------
    pandas dataframe
        initial dataframe in which all lines ARE NOT in the range of each entry of the ROD dictionary.
    """
    if not ROD:
        return df
    for key, value in ROD.items():
        # Rappel : key est par ex "Vx" ; value est {"size":10, "position":0}
        if "range" in value:
            minimum = np.min(value["range"])
            maximum = np.max(value["range"])
        elif "minimum" in value and "maximum" in value:
            minimum = value["minimum"]
            maximum = value["maximum"]
        elif "min" in value and "max" in value:
            minimum = value["min"]
            maximum = value["max"]
        elif type(value) == list or type(value) == tuple:
            minimum = min(value)
            maximum = max(value)
        else:
            print("[WARNING] The ROD function was not recognized. Please fix me.")
            return df

        if key in df:
            df = df[~((df[key] >= minimum) & (df[key] < maximum))]
        else:
            print(f"[WARNING] The key {key} of the ROI is not in the other dataframe.")
    return df


def apply_ROI(atoms, ROI):
    """
    This function returns athe atoms dataframe such that all elements are within the ROI dictionary range. If the ROI is an empty dictionary, this function returns the initial dataframe.

    Parameters
    ----------
    atoms : pandas dataframe
        DataFrame with all atoms.
    ROI : dic
        dictionary for which every entry (for exemple 'T') matches a column of the atoms dataframe. The function returns a dictionary with the same number of columns as atoms but for which every line is in a range required by the ROI dictionary, i.e. a maximum and a minimum value. Since a recent update, each entry of the dictionary can be a tuple or a list with two number or a dictionary with entries "min/max" or 'range'.

    Returns
    -------
    pandas dataframe
        initial dataframe in which all lines ARE in the range of each entry of the ROI dictionary.
    """
    if ROI:
        for key, value in ROI.items():
            # Rappel : key est par ex "Vx" ; value est {"size":10, "position":0}
            if "range" in value:
                minimum = np.min(value["range"])
                maximum = np.max(value["range"])
            elif "minimum" in value and "maximum" in value:
                minimum = value["minimum"]
                maximum = value["maximum"]
            elif "min" in value and "max" in value:
                minimum = value["min"]
                maximum = value["max"]
            elif type(value) == list or type(value) == tuple:
                minimum = min(value)
                maximum = max(value)
            else:
                print(
                    "[WARNING] The ROI format was not recognized. Please read the apply_ROI documentation. We expect a dictionary with all entrien"
                )
                return atoms
            if key in atoms:
                atoms = atoms[((atoms[key] <= maximum) & (atoms[key] > minimum))]
            else:
                print(
                    f"[WARNING] The key {key} of the ROI is not in the other dataframe."
                )
    return atoms


def check_BEC_fit(
    folder,
    histogramm_width=0.01,
    ROI_for_fit={
        "T": {"min": 305.5, "max": 309.7},
        "X": {"min": -35, "max": -7},
        "Y": {"min": -20, "max": 20},
    },
    width_saturation=0,
):
    atom_files = select_atoms_in_folder(folder)
    X, Y, T = load_XYTTraw(atom_files[random.randint(0, len(atom_files))])
    fit_BEC_arrival_time(
        X,
        Y,
        T,
        show_fit=True,
        histogramm_width=histogramm_width,
        ROI_for_fit=ROI_for_fit,
        width_saturation=width_saturation,
    )


def fit_BEC_arrival_time(
    X,
    Y,
    T,
    ROI_for_fit={
        "T": {"min": 305.5, "max": 309.7},
        "X": {"min": -35, "max": -7},
        "Y": {"min": -35, "max": 35},
    },
    histogramm_width=0.01,
    width_saturation=0,
    show_fit=False,
):
    """
    This functions fit BEC arrival times. It generates a dictionary named ans in which we store some properties of the arrival times of our BEC.
    """
    ans = {"Number of Atoms": len(X)}
    data = pd.DataFrame({"X": X, "Y": Y, "T": T})
    data = apply_ROI(data, ROI_for_fit)

    X = data["X"].to_numpy()
    Y = data["Y"].to_numpy()
    T = data["T"].to_numpy()
    ans["Number of Atoms in ROIfit"] = len(T)
    ans["BEC Std Arrival Time"] = np.std(T)

    #### FIT IN TIME

    bin_heights, bin_borders = np.histogram(
        T,
        bins=np.arange(
            ROI_for_fit["T"]["min"], ROI_for_fit["T"]["max"], histogramm_width
        ),
    )
    bin_centers = np.array(bin_borders[:-1] + np.diff(bin_borders) / 2)
    # find the position of the max
    max_index = np.argmax(bin_heights)
    mean = bin_centers[max_index]
    ans["BEC Arrival Time with max"] = mean
    sigma = np.mean(bin_heights * (bin_centers - mean) ** 2)
    sigma = 0.1
    p0 = [mean, np.max(bin_heights), sigma, 0]
    bin_heights = list(bin_heights)
    bin_centers = list(bin_centers)
    # ci -dessous, on supprime un certain nombre de points pour ne pas prendre en compte la saturation du mcp.
    n_hole = int(width_saturation / histogramm_width)
    failed_status = False
    for i in range(n_hole):
        bin_centers.pop(max_index + 1)
        bin_heights.pop(max_index + 1)
    try:
        popt, pcov = curve_fit(gaussian_function, bin_centers, bin_heights, p0=p0)
        # perr = np.sqrt(np.diag(pcov))
    except:
        failed_status = True
        popt = p0
        # perr = [np.nan, np.nan, np.nan]
    ans["BEC Arrival Time"] = popt[0]
    ans["BEC fitted Std Arrival Time"] = popt[2]
    ans["BEC Arrival Time with fit"] = popt[0]
    ### FIT IN X
    bin_heightsX, bin_bordersX = np.histogram(
        X,
        bins=np.arange(ROI_for_fit["X"]["min"], ROI_for_fit["X"]["max"]),
    )
    bin_centersX = np.array(bin_bordersX[:-1] + np.diff(bin_bordersX) / 2)
    # find the position of the max
    max_indexX = np.argmax(bin_heightsX)
    arr_time_maximumX = bin_centersX[max_indexX]
    meanX = bin_centersX[max_indexX]
    offset = np.min(bin_heightsX)
    sigmaX = np.mean((bin_heightsX - offset) * (bin_centersX - meanX) ** 2)
    sigmaX = 7
    p0X = [np.max(meanX), np.max(bin_heightsX), np.abs(sigmaX), 0]
    bin_heightsX = list(bin_heightsX)
    bin_centersX = list(bin_centersX)
    try:
        poptX, pcovX = curve_fit(gaussian_function, bin_centersX, bin_heightsX, p0=p0X)
        # perr = np.sqrt(np.diag(pcov))
    except:
        failed_status = True
        poptX = p0X
    ans["BEC Center X"] = poptX[0]
    ### FIT in Y
    bin_heightsY, bin_bordersY = np.histogram(
        Y,
        bins=np.arange(ROI_for_fit["Y"]["min"], ROI_for_fit["Y"]["max"]),
    )
    bin_centersY = np.array(bin_bordersY[:-1] + np.diff(bin_bordersY) / 2)
    # find the position of the max
    max_indexY = np.argmax(bin_heightsY)
    arr_time_maximumY = bin_centersY[max_indexY]
    meanY = bin_centersY[max_indexY]
    offset = np.min(bin_heightsY)
    sigmaY = np.mean((bin_heightsY - offset) * (bin_centersY - meanY) ** 2)
    sigmaY = 7
    p0Y = [np.max(meanY), np.max(bin_heightsY), np.abs(sigmaY), offset]
    bin_heightsY = list(bin_heightsY)
    bin_centersY = list(bin_centersY)
    try:
        poptY, pcovY = curve_fit(gaussian_function, bin_centersY, bin_heightsY, p0=p0Y)
        # perr = np.sqrt(np.diag(pcov))
    except:
        failed_status = True
        poptY = p0Y
    ans["BEC Center Y"] = poptY[0]
    if show_fit:
        fig, axes = plt.subplots(figsize=(10, 3), ncols=3)
        axes[0].plot(bin_centers, bin_heights, "*", label="data")
        print(" ##########  FIT   ##########")
        print("p0 : [mean, amplitude, standard_deviation, offset]")
        print("Fit in T :")
        print(f"p0 : {p0}")
        print(f"popt : {popt}")
        print("=" * 20)
        axes[0].plot(
            bin_centers, gaussian_function(bin_centers, *popt), "-", label="fit"
        )
        axes[0].plot(
            bin_centers, gaussian_function(bin_centers, *p0), "--", label="guess"
        )
        print("Fit in X :")
        print(f"p0 : {p0X}")
        print(f"popt : {poptX}")
        print("=" * 20)
        axes[1].plot(bin_centersX, bin_heightsX, "*", label="data")
        axes[1].plot(
            bin_centersX, gaussian_function(bin_centersX, *poptX), "-", label="fit"
        )
        axes[1].plot(
            bin_centersX, gaussian_function(bin_centersX, *p0X), "--", label="guess"
        )
        print("Fit in X Y")
        print(f"p0 : {p0Y}")
        print(f"popt : {poptY}")
        print("=" * 20)
        axes[2].plot(bin_centersY, bin_heightsY, "*", label="data")
        axes[2].plot(
            bin_centersY, gaussian_function(bin_centersY, *poptY), "-", label="fit"
        )
        axes[2].plot(
            bin_centersY, gaussian_function(bin_centersY, *p0Y), "--", label="guess"
        )
        plt.show()

    return ans, failed_status


def obtain_arrival_times(atom_files, **kwargs):
    """Pour chaque cycle du dataframe atoms, on fait un histogramme des temps d'arrivée des atomes puis on considère que le max de cet histogramme correspond au temps d'arrivé du BEC. Le paramètre histogramm_width permet d'ajuster la largeur des pics de l'histogramme (ms)
    WARNING : cette fonction est longue et non optimisée.

    Parameters
    ----------
    atoms : list of string path
        liste des cycles avec les chemins allant au .atoms.
    **kwargs : any supplementary keyword argument is given to the fit_BEC_arrival_time function.

    Return
    ------
    dataframe containing BEC arrival time & # of atoms
    """
    print("Starting to gather arrival time of BECs")
    list_of_df = []
    for i, path in enumerate(tqdm(atom_files)):
        X, Y, T = load_XYTTraw(path)
        seq, cycle = return_cycle_from_path(path)
        ans, failed_status = fit_BEC_arrival_time(X, Y, T, show_fit=False, **kwargs)
        ans["Cycle"] = cycle
        list_of_df.append(pd.DataFrame(ans, index=[i]))
        df_arrival_time = pd.concat(list_of_df)
        list_of_df = [df_arrival_time]

    return df_arrival_time


def export_data_set_to_pickle(
    folder,
    ROI,
    ROD={"T": [10, 20]},
    find_arrival_times=False,
    n_max_cycles=1e8,
    histogramm_width=0.01,
    ROI_for_fit={"T": {"min": 306.2, "max": 309.7}},
    width_saturation=0,
    supplementary_rois=[],
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
    atoms_in_ROI = apply_ROD(atoms_in_ROI, ROD)
    filename_dataset = os.path.join(folder, "dataset.pkl")
    atoms_in_ROI.to_pickle(filename_dataset)
    df_parameters = gather_saved_sequence_parameters(folder)
    message = "Warning: since a new version of export_dataset_to_pickle, we fit the X and Y position of the BEC and a ROI should be given in the ROI_for_fit dictionary."
    if find_arrival_times:
        if "X" not in ROI_for_fit.keys():
            ROI_for_fit["X"] = {"min": -35, "max": -7}
            print(message)
        if "Y" not in ROI_for_fit.keys():
            ROI_for_fit["Y"] = {"min": -35, "max": 35}
            print(message)
        df_arrival_times = obtain_arrival_times(
            atom_files,
            histogramm_width=histogramm_width,
            ROI_for_fit=ROI_for_fit,
            width_saturation=width_saturation,
        )

        df_arrival_times = pd.merge(
            df_arrival_times,
            atoms_in_ROI.groupby("Cycle")
            .count()["T"]
            .rename("Number of Atoms in ROI")
            .reset_index(),
            on="Cycle",
        )
        for i, roi in enumerate(supplementary_rois):
            at = apply_ROI(atoms, roi)
            df_arrival_times = pd.merge(
                df_arrival_times,
                at.groupby("Cycle")
                .count()["T"]
                .rename("Number of Atoms in supplementary ROI{}".format(i))
                .reset_index(),
                on="Cycle",
            )
        df_arrival_times = pd.merge(df_arrival_times, df_parameters, on="Cycle")
        filename = os.path.join(folder, "arrival_times.pkl")
        df_arrival_times.to_pickle(filename)

    return filename_dataset


def gather_HAL_fits(folder):
    folder = Path(folder + "/.HAL_fits")
    ext = [".json"]
    path_list = sorted(
        [
            path.as_posix()
            for path in filter(lambda path: path.suffix in ext, folder.glob("*"))
        ]
    )  # string object
    dataframe = pd.DataFrame()
    for filename in path_list:
        sequence_parameters_filename = filename.replace(".HAL_fits/", "")
        seq, cycle = return_cycle_from_path(sequence_parameters_filename)
        column_names = ["Sequence", "Cycle"]
        column_values = [seq, cycle]
        if os.path.exists(filename):
            f = open(filename)
            data = json.load(f)
            data = data["collection"]
            for ROI_name, ROI in data.items():
                for key, values in ROI["result"].items():
                    if key == "values":
                        for element in values:
                            column_names.append(
                                "{} | {} ({})".format(
                                    ROI_name, element["name"], element["unit"]
                                )
                            )
                            column_values.append(element["value"])
        new_df = pd.DataFrame([column_values], columns=column_names)
        dataframe = pd.concat([dataframe, new_df])
    dataframe.reset_index(drop=True, inplace=True)
    filename_parameters = os.path.join(folder, "hal_fits.pkl")
    dataframe.to_pickle(filename_parameters)
    return dataframe


def gather_saved_sequence_parameters(folder):
    """Cette fonction récupère tous les paramètres sauvegardés pour chaque cycle de
    la séquence et les mets dans un dataframe avec le numéro du cycle.

    Parameters
    ----------
    folder : path like
        chemin vers le dossier contenant tous les .atoms
    """
    folder = Path(folder)
    ext = [".json"]
    path_list = sorted(
        [
            path.as_posix()
            for path in filter(lambda path: path.suffix in ext, folder.glob("*"))
        ]
    )  # string object
    dataframe = pd.DataFrame()
    for filename in path_list:
        # filename = atom_name.replace(".atoms", ".json")
        seq, cycle = return_cycle_from_path(filename)
        column_names = ["Sequence", "Cycle"]
        column_values = [seq, cycle]
        if os.path.exists(filename):
            f = open(filename)
            data = json.load(f)

            for element in data:
                if element["name"] not in [
                    "sequence",
                    "cycle",
                    "cycle_id",
                    "Sequence",
                    "Cycle",
                ]:
                    column_names.append(
                        "{} ({})".format(element["name"], element["unit"])
                    )
                    column_values.append(element["value"])
        new_df = pd.DataFrame([column_values], columns=column_names)
        dataframe = pd.concat([dataframe, new_df])
    dataframe.reset_index(drop=True, inplace=True)
    filename_parameters = os.path.join(folder, "parameters.pkl")
    dataframe.to_pickle(filename_parameters)
    return dataframe


if __name__ == "__main__":
    # folder = Path("/mnt/manip_E/2022/11/18/083")
    # export_data_set_to_pickle(
    #     folder,
    #     ROI={"T": {"min": 306.2, "max": 309.7}},
    #     find_arrival_times=True,
    #     n_max_cycles=3,
    # )
    folder = "/mnt/manip_E/2023/06/16/060"
    gather_saved_sequence_parameters(folder)
    # folder = "/media/victor/8482-1010/gus_data/2023/05/15/053"
    # folder = "/home/victor/gus_data/2023/05/15/053"
    # check_BEC_fit(
    #     folder,
    #     width_saturation=0.3,
    #     histogramm_width=0.01,
    #     ROI_for_fit={
    #         "T": {"min": 307, "max": 309.7},
    #         "X": {"min": -35, "max": -7},
    #         "Y": {"min": -35, "max": 35},
    #     },
    # )

    # X, Y, T = load_XYTTraw(folder + "/041_225.atoms")
    # print(
    #     fit_BEC_arrival_time(
    #         T,
    #         show_fit=True,
    #         width_saturation=0.5,
    #         ROI_for_fit={"T": {"min": 307, "max": 309.7}},
    #     )
    # )

    # atom_files = select_atoms_in_folder(folder)
    # df_parameters = gather_saved_sequence_parameters(folder)
    # df_arrival_times = obtain_arrival_times(
    #     atom_files,
    #     histogramm_width=0.01,
    #     width_saturation=0.2,
    # )
    # df_arrival_times = pd.merge(df_arrival_times, df_parameters, on="Cycle")
    # filename = os.path.join(folder, "arrival_times.pkl")
    # df_arrival_times.to_pickle(filename)
    # folder = "/media/victor/8482-1010/gus_data/2023/05/15/041"
    # # X, Y, T = load_XYTTraw(folder + "/041_225.atoms")
    # # print(
    # #     fit_BEC_arrival_time(
    # #         T,
    # #         show_fit=True,
    # #         width_saturation=0.5,
    # #         ROI_for_fit={"T": {"min": 307, "max": 309.7}},
    # #     )
    # # )

    # atom_files = select_atoms_in_folder(folder)
    # df_parameters = gather_saved_sequence_parameters(folder)
    # df_arrival_times = obtain_arrival_times(
    #     atom_files,
    #     histogramm_width=0.01,
    #     width_saturation=0.2,
    # )
    # df_arrival_times = pd.merge(df_arrival_times, df_parameters, on="Cycle")
    # filename = os.path.join(folder, "arrival_times.pkl")
    # df_arrival_times.to_pickle(filename)
