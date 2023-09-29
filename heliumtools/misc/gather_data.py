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
    print(
        "[WARNING] : This function should be deleted in the next version of Heliumtools. Please use apply_ROI(df -> pd.DataFrame, ROI -> dictionary) instead."
    )
    return apply_ROI(df, roi)


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
    time_resolution = 2.75e-10  # old
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


def load_times_of_atoms(path):
    """returns X1, X2, Y1, Y2 of a reconstructed atom.

    Parameters
    ----------
    path : str or PathLike object
        path to the .atoms file
    """
    atoms = np.fromfile(path, dtype="uint64")
    events_list = atoms.reshape(int(len(atoms) / 4), 4).T
    return events_list[0], events_list[1], events_list[2], events_list[3]


def load_offset_of_atoms(path):
    path = str(path)
    path = path.replace(".atoms", ".offsets")
    offsets = np.fromfile(path, dtype="int64")
    return offsets


def load_raw_time(times_path):
    # if the .atoms is from a rereconstruction, there are no .atoms files hence we must recover them using the other algorithm.
    if ~os.path.exists(times_path):
        config_path = times_path.split(".times")[0] + ".stats"
        myfile = open(config_path, "r")
        lines = myfile.readlines()

        for index, line in enumerate(lines):
            content = line.strip()
            if "Origin file : " in content:
                times_path = (
                    content.replace("Origin file : ", "")
                    + ".times"
                    + times_path.split(".times")[1]
                )
                break

        myfile.close()
    times = np.fromfile(times_path, dtype="uint64")
    time_resolution = 1.2e-10
    times = times * time_resolution * 1e3
    return times


def data_filter(data, bec_arrival_times, filters):
    selec_bec_arrival_times = apply_ROI(bec_arrival_times, filters)
    selected_data = data[data["Cycle"].isin(selec_bec_arrival_times["Cycle"])]
    return selected_data, selec_bec_arrival_times


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
        (minimum, maximum) = get_roi_min_max(ROD, key)

        if key in df:
            df = df[~((df[key] >= minimum) & (df[key] < maximum))]
        else:
            print(f"[WARNING] The key {key} of the ROI is not in the other dataframe.")
    return df


def apply_ROI(atoms, ROI):
    """
    This function returns the atoms dataframe such that all elements are within the ROI dictionary range. If the ROI is an empty dictionary, this function returns the initial dataframe.

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
            (minimum, maximum) = get_roi_min_max(ROI, key)
            if key in atoms:
                atoms = atoms[((atoms[key] <= maximum) & (atoms[key] > minimum))]
            else:
                print(
                    f"[WARNING] The key {key} of the ROI is not in the other dataframe."
                )
    return atoms


def get_roi_min_max(roi, axis):
    """This function returns the maximum and minimum value of a roi type dictionary for a given entry axis.

    Parameters
    ----------
    roi : dictionary
        ROI type dictionary
    axis : str
        key of the dictionary from which we we want the minimum and maximum value

    Returns
    -------
    tuples
        maximum and minimum value of the key axis of the roi.
    """
    if axis not in roi:
        print(f"[WARNING] : the axis {axis} is not in the ROI.")
        return (-np.inf, np.inf)
    value = roi[axis]
    if "range" in value:
        minimum = np.min(value["range"])
        maximum = np.max(value["range"])
    elif "minimum" in value and "maximum" in value:
        minimum = value["minimum"]
        maximum = value["maximum"]
    elif "min" in value and "max" in value:
        minimum = value["min"]
        maximum = value["max"]
    elif "position" in value and "size" in value:
        minimum = value["position"] - 0.5 * value["size"]
        maximum = value["position"] + 0.5 * value["size"]
    elif type(value) == list or type(value) == tuple:
        minimum = min(value)
        maximum = max(value)
    else:
        print(
            "[WARNING] The ROI format was not recognized. Please read the apply_ROI documentation. We expect a dictionary with all values being either a dictionary or a list. "
        )
    return (minimum, maximum)


def get_roi_size(roi, axis):
    """Returns the size of a ROI like dictionary along a given axis."""
    (minimum, maximum) = get_roi_min_max(roi, axis)
    return maximum - minimum


def check_roi_for_fit(roi):
    """Check if the roi use for the fit of the BEC arrival time has enough entries."""
    if "T" not in roi:
        print(
            f"[WARNING] There is no entry T in the ROI for fit. Adding default value [305.5, 309.7]."
        )
        roi["T"] = [305.5, 309.7]
    if "X" not in roi:
        print(
            f"[WARNING] There is no entry X in the ROI for fit. Adding default value [-35,-7]."
        )
        roi["X"] = [-35, -7]
    if "Y" not in roi:
        print(
            f"[WARNING] There is no entry Y in the ROI for fit. Adding default value [-35, 35]."
        )
        roi["Y"] = [-35, 35]
    return roi


def fit_BEC_arrival_time(
    data,
    filename,
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
    ans = {"Number of Atoms": len(data)}
    ROI_for_fit = check_roi_for_fit(ROI_for_fit)
    data = apply_ROI(data, ROI_for_fit)
    X = data["X"].to_numpy()
    Y = data["Y"].to_numpy()
    T = data["T"].to_numpy()
    ans["Number of Atoms in ROIfit"] = len(data)
    ans["BEC Std Arrival Time"] = np.std(data["T"])
    if show_fit:
        fig, axes = plt.subplots(figsize=(3.3 * 4, 3 * 2), ncols=4, nrows=2)
        axes = axes.flatten()
        print(" ##########  FIT   ##########")
        print("p0 : [mean, amplitude, standard_deviation, offset]")

    #### FIT IN TIME
    if True:
        mini, maxi = get_roi_min_max(ROI_for_fit, "T")
        bin_heights, bin_borders = np.histogram(
            T,
            bins=np.arange(mini, maxi, histogramm_width),
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
            if max_index + 2 < len(bin_centers):
                bin_centers.pop(max_index + 1)
                bin_heights.pop(max_index + 1)
        try:
            popt, pcov = curve_fit(gaussian_function, bin_centers, bin_heights, p0=p0)
            # perr = np.sqrt(np.diag(pcov))
        except:
            failed_status = True
            popt = p0
        ans["BEC Arrival Time"] = popt[0]
        ans["BEC fitted Std Arrival Time"] = popt[2]
        ans["BEC Arrival Time with fit"] = popt[0]
        if show_fit:
            ax = axes[0]
            ax.plot(bin_centers, bin_heights, "o", alpha=0.8, label="data")
            print("Fit in T :")
            print(f"p0 : {p0}")
            print(f"popt : {popt}")
            print("=" * 20)
            ax.plot(
                bin_centers, gaussian_function(bin_centers, *popt), "-", label="fit"
            )
            ax.plot(
                bin_centers, gaussian_function(bin_centers, *p0), "--", label="guess"
            )
            ax.set_title("Mean : {:.3f} ms".format(popt[0]))
            ax.set_xlabel("Arrival time of reconstructed atoms (ms)")
    ##### FIT in X and Y
    for i, XY in enumerate(["X", "Y"]):
        (bin_mini, bin_maxi) = get_roi_min_max(ROI_for_fit, XY)
        bin_heightsXY, bin_bordersXY = np.histogram(
            data[XY].to_numpy(),
            bins=np.arange(bin_mini, bin_maxi),
        )
        bin_centersXY = np.array(bin_bordersXY[:-1] + np.diff(bin_bordersXY) / 2)
        # find the position of the max
        max_indexXY = np.argmax(bin_heightsXY)
        arr_time_maximumXY = bin_centersXY[max_indexXY]
        mean = bin_centersXY[max_indexXY]
        offset = np.min(bin_heightsXY)
        # sigma = np.mean((bin_heightsXY - offset) * (bin_centersXY - mean) ** 2)
        sigma = 4
        p0XY = [mean, np.max(bin_heightsXY) - offset, np.abs(sigma), offset]
        bin_heightsXY = list(bin_heightsXY)
        bin_centersXY = list(bin_centersXY)
        try:
            poptXY, pcovXY = curve_fit(
                gaussian_function, bin_centersXY, bin_heightsXY, p0=p0XY
            )
        except:
            failed_status = True
            poptXY = p0XY
        ans["BEC Center " + XY] = poptXY[0]
        if show_fit:
            print("Fit in " + XY + " :")
            print(f"p0 : {p0XY}")
            print(f"popt : {poptXY}")
            print("=" * 20)
            ax = axes[1 + i]
            ax.plot(bin_centersXY, bin_heightsXY, "*", alpha=0.7, label="data")
            ax.plot(
                bin_centersXY,
                gaussian_function(bin_centersXY, *poptXY),
                "-",
                label="fit",
            )
            ax.plot(
                bin_centersXY,
                gaussian_function(bin_centersXY, *p0XY),
                "--",
                label="guess",
            )
            ax.set_title("Mean : {:.3f} mm".format(poptXY[0]))
            ax.set_xlabel(XY + " position of reconstructed atoms (mm)")
    ##### FIT of each channel X1, x2, y1 and y2
    if filename:
        ans["Mean Arrival Time (fit .times)"] = 0
        for index, xj in enumerate(["x1", "x2", "y1", "y2"]):
            try:
                times = load_raw_time(filename.replace(".atoms", ".times" + xj))
            except:
                break
            mini, maxi = get_roi_min_max(ROI_for_fit, "T")
            times = times[np.where(times > mini)]
            times = times[np.where(times < maxi)]
            ans["BEC Std of " + xj] = np.std(times)

            bin_heights, bin_borders = np.histogram(
                times,
                bins=np.arange(mini, maxi, histogramm_width),
            )
            bin_centers = np.array(bin_borders[:-1] + np.diff(bin_borders) / 2)
            # find the position of the max
            max_index = np.argmax(bin_heights)
            mean = bin_centers[max_index]
            ans[xj + " Arrival Time (maximum)"] = mean
            sigma = np.mean(bin_heights * (bin_centers - mean) ** 2)
            sigma = 0.1
            p0 = [mean, np.max(bin_heights), sigma, 0]
            bin_heights = list(bin_heights)
            bin_centers = list(bin_centers)
            # ci -dessous, on supprime un certain nombre de points pour ne pas prendre en compte la saturation du mcp.
            n_hole = int(width_saturation / histogramm_width)
            failed_status = False
            for i in range(n_hole):
                if max_index + 2 < len(bin_centers):
                    bin_centers.pop(max_index + 1)
                    bin_heights.pop(max_index + 1)
            try:
                popt, pcov = curve_fit(
                    gaussian_function, bin_centers, bin_heights, p0=p0
                )
                # perr = np.sqrt(np.diag(pcov))
            except:
                failed_status = True
                popt = p0
            ans[xj + " Arrival Time (fit)"] = popt[0]
            ans["Mean Arrival Time (fit .times)"] += popt[0]
            if show_fit:
                ax = axes[3 + index]
                ax.plot(bin_centers, bin_heights, "o", alpha=0.8, label="data")
                print("Fit in T :")
                print(f"p0 : {p0}")
                print(f"popt : {popt}")
                print("=" * 20)
                ax.plot(
                    bin_centers, gaussian_function(bin_centers, *popt), "-", label="fit"
                )
                ax.plot(
                    bin_centers,
                    gaussian_function(bin_centers, *p0),
                    "--",
                    label="guess",
                )
                ax.set_title("Mean : {:.3f} ms".format(popt[0]))
                ax.set_xlabel("Unreconstructed signal " + xj)
        ans["Mean Arrival Time (fit .times)"] = (
            ans["Mean Arrival Time (fit .times)"] / 4
        )
    if show_fit:
        for ax in axes:
            ax.legend(
                loc=0,
            )
        plt.tight_layout()
        plt.show()
    return ans, failed_status


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
    selected_file = atom_files[random.randint(0, len(atom_files))]
    X, Y, T = load_XYTTraw(selected_file)
    data = pd.DataFrame({"X": X, "Y": Y, "T": T})
    ans = fit_BEC_arrival_time(
        data,
        selected_file,
        show_fit=True,
        histogramm_width=histogramm_width,
        ROI_for_fit=ROI_for_fit,
        width_saturation=width_saturation,
    )

    print(ans)


def fit_BEC_arrival_time_old(
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
    ROI_for_fit = check_roi_for_fit(ROI_for_fit)
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
        if max_index + 2 < len(bin_centers):
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
        print("Fit in " + xj + " :")
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
    selected_files = select_atoms_in_folder(folder)
    ROI_for_fit = check_roi_for_fit(ROI_for_fit)
    N_files = min([len(selected_files), n_max_cycles])
    # on va remplir les listes suivantes avec l'ensemble des atomes.
    Xc, Yc, Tc, Cyclec = [], [], [], []
    selected_files = selected_files[0:N_files]
    print("Starting to gather atoms")
    df_atoms = pd.DataFrame()
    df_parameters = pd.DataFrame()
    df_arrival_times = pd.DataFrame()
    for i, filename in tqdm(enumerate(selected_files)):
        X, Y, T = load_XYTTraw(filename)
        seq, cycle = return_cycle_from_path(filename)
        raw_data = pd.DataFrame({"X": X, "Y": Y, "T": T})
        raw_data["Cycle"] = cycle
        atoms_in_ROI = apply_ROI(raw_data, ROI)
        atoms_in_ROI = apply_ROD(atoms_in_ROI, ROD)
        df_atoms = pd.concat([df_atoms, atoms_in_ROI])

        ## On charge maintenant les parametres de séquence
        # filename = atom_name.replace(".atoms", ".json")
        seq_par = get_sequence_parameters(filename.replace(".atoms", ".json"))
        seq_par["Sequence"] = seq
        seq_par["Cycle"] = cycle
        new_df = pd.DataFrame(seq_par, index=[i])
        df_parameters = pd.concat([df_parameters, new_df])

        ## Fits des temps d'arrivée
        if find_arrival_times:
            arrival, failed_status = fit_BEC_arrival_time(
                raw_data,
                filename,
                show_fit=False,
                histogramm_width=histogramm_width,
                ROI_for_fit=ROI_for_fit,
                width_saturation=width_saturation,
            )
            arrival["Cycle"] = cycle
            arrival["Number of Atoms in ROI"] = len(atoms_in_ROI)
            arrival["dT of atoms in ROI"] = np.std(atoms_in_ROI["T"])
            for i, roi in enumerate(supplementary_rois):
                at = apply_ROI(raw_data, roi)
                arrival["Number of Atoms in supplementary ROI{}"] = len(at)

            df_arrival_times = pd.concat(
                [df_arrival_times, pd.DataFrame(arrival, index=[i])]
            )

    filename_dataset = os.path.join(folder, "dataset.pkl")
    df_atoms.to_pickle(filename_dataset)
    filename_parameters = os.path.join(folder, "parameters.pkl")
    df_parameters.to_pickle(filename_parameters)
    if find_arrival_times:
        filename = os.path.join(folder, "arrival_times.pkl")
        df_arrival_times = df_arrival_times.merge(df_parameters, on="Cycle")
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


def get_sequence_parameters(filename):
    seq_par = {}
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
                seq_par["{} ({})".format(element["name"], element["unit"])] = element[
                    "value"
                ]

    return seq_par


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
    for i, filename in enumerate(path_list):
        # filename = atom_name.replace(".atoms", ".json")
        seq, cycle = return_cycle_from_path(filename)
        seq_par = get_sequence_parameters(filename)
        seq_par["Sequence"] = seq
        seq_par["Cycle"] = cycle
        new_df = pd.DataFrame(seq_par, index=[i])
        dataframe = pd.concat([dataframe, new_df])
    print(dataframe)
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

    if True:  # check du fit du BEC
        folder = "/mnt/manip_E/2023/07/12/030"
        folder = "/mnt/manip_E/2023/08/15/006"
        check_BEC_fit(
            folder,
            width_saturation=0.3,
            histogramm_width=0.01,
            ROI_for_fit={
                "T": {"min": 305, "max": 310},
                "X": {"min": -35, "max": -7},
            },
        )
    if False:
        folder = "/mnt/manip_E/2023/08/16/006"
        df_parameters = gather_saved_sequence_parameters(folder)
    if False:
        folder = "/mnt/manip_E/2023/08/16/006"
        export_data_set_to_pickle(
            folder, ROI={"T": [312, 322]}, find_arrival_times=True
        )
    # folder = "/media/victor/8482-1010/gus_data/2023/05/15/053"
    # folder = "/home/victor/gus_data/2023/05/15/053"

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
