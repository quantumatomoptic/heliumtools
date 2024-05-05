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
import glob, os, re, json, random, platform
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
from tqdm import tqdm
from heliumtools.misc import logger
import warnings
from scipy.optimize import OptimizeWarning
from flatten_dict import flatten, unflatten, reducers, splitters
import pytz
from datetime import datetime
from heliumtools.tools import (
    apply_ROD,
    apply_ROI,
    data_filter,
    get_roi_min_max,
    get_roi_size,
    get_roi_center,
)

log = logger.getLogger(__name__)


def gaussian_function(x, mean, amplitude, standard_deviation, offset):
    return offset + amplitude * np.exp(-((x - mean) ** 2) / (2 * standard_deviation**2))


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


def return_cycle_from_path(path, show_error=True):
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
        elif "." in path:
            a, b = path.split(".")
            json_path = a + ".json"
        else:
            json_path = path + ".json"
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
        msg = f"Failed to find the cycle of path {path} --> I try to guess it from the filename."
        if show_error:
            log.warning(msg)
    ## If the cycle was not detected, we interpret it
    ## with the filename
    try:
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
    except:
        log.error(f"No cycle was identified in cycle {path}")
        return 0, 0


def load_atoms(folder, n_max_cycles=1e8):
    """Load all .atoms files

    Parameters
    ----------
    folder : path like
        chemin vers le dossier contenant tous les .atoms
    """
    selected_files = select_atoms_in_folder(Path(folder))
    N_files = min([len(selected_files), n_max_cycles])
    Xc, Yc, Tc, Cyclec = [], [], [], []
    selected_files = selected_files[0:N_files]
    print(f"Starting to gather atoms from {folder}")
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
    """Return the X, Y, T positions of a cycle in a sequence

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

    Parameters
    ----------
    data : pandas DataFrame
        DataFrame with 3 columns X, Y, T, the data to fit.
    filename : str
        Name of the file from which the data comes. This name will be used to retrieve the corresponding .times file for this cycle.
    ROI_for_fit : dict
        Dictionary with the ROI for fitting. This ROI is applied to all axes (X, Y, T, and .times).
    histogramm_width : float
        Width of the histogram bins for fitting according to T and .times, in ms.
    width_saturation : float
        Saturation width during which there is no signal due to TDC saturation. The histogram points between tmax, the time at which the signal is maximal, and tmax + dt are removed and not taken into account in the fit.
    show_fit : bool
        Whether to display the fit on a figure (do not set to True if performing multiple fits). Consider using the check_BEC_fit function if you want to verify your fits.

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
            ax.axvspan(
                bin_centers[max_index],
                bin_centers[max_index] + n_hole * histogramm_width,
                alpha=0.2,
                color="red",
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
                ax.axvspan(
                    bin_centers[max_index],
                    bin_centers[max_index] + n_hole * histogramm_width,
                    alpha=0.2,
                    color="red",
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
    """fit function that checks if the BEC fit is correct or not.

    Parameters
    ----------
    folder : str or path like
        folder in which one wats to extract a cycle.
    histogramm_width : float
    Width of the histogram bins for fitting according to T and .times, in ms.
    ROI_for_fit : dict, optional
        Dictionary with the ROI for fitting. This ROI is applied to all axes (X, Y, T, and .times)., by default { "T": {"min": 305.5, "max": 309.7}, "X": {"min": -35, "max": -7}, "Y": {"min": -20, "max": 20}, }
    width_saturation : float
        Saturation width during which there is no signal due to TDC saturation. The histogram points between tmax, the time at which the signal is maximal, and tmax + dt are removed and not taken into account in the fit.


    """
    atom_files = select_atoms_in_folder(folder)
    selected_file = atom_files[random.randrange(0, len(atom_files))]
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


def obtain_arrival_times(atom_files, **kwargs):
    """For each cycle of the 'atoms' dataframe, we create a histogram of atom arrival times, then we consider the maximum of this histogram as the arrival time of the BEC. The 'histogram_width' parameter adjusts the width of the histogram peaks (in ms).
    WARNING: This function is lengthy and not optimized

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


def get_sequence_id_from_seq_dir(seq_dir: str):
    """returns the sequence id  string and integer from the sequence directory folder

    Parameters
    ----------
    seq_dir : str
        sequence directory string

    Returns
    -------
    str
        sequence id string
    int
        sequence id number
    """
    try:
        seq_dir = str(seq_dir)
        seq_id = seq_dir.replace("/", "-").replace("\\", "-").replace(os.path.sep, "-")
        match = re.search(r"\b\d{4}-\d{2}-\d{2}-\d{3}\b", seq_id)
        if match:
            seq_id = match.group(0)  # 2024-08-02-005 for example
            seq_num = int(seq_id.replace("-", "")[2:])  # 240802005 in this case

            return seq_id, seq_num
        else:
            msg = f"The sequence ID was not recognized from {seq_id}"
            log.error(msg)
            return seq_dir, 0
    except:
        msg = f"The sequence ID was not recognized from {seq_dir}"
        log.warning(msg)
        return seq_dir, 0


def export_data_set_to_pickle(
    folder,
    ROI,
    ROD={"T": [10, 20]},
    find_arrival_times=False,
    n_max_cycles=1e8,
    histogramm_width=0.01,
    ROI_for_fit={"T": {"min": 306.2, "max": 309.7}, "X": [-30, 0], "Y": [-15, 15]},
    width_saturation=0,
    supplementary_rois=[],
    metadata=[],
):
    """
    Export the dataset folder as a pickle. This function generates three files in the folder. Thus, the dataframe folder/dataset.pkl contains all the atoms from the folder within the ROI. This function also creates folder/parameters.pkl containing all the sequence parameters, as well as folder/arrival_times.pkl containing the fitted arrival times (plus parameters) of the sequence.

    Parameters
    ----------
    folder : path-like
        Path to the folder containing all the .atoms files.
    ROI : dict
        Region of Interest: we will select all atoms that are within this ROI. It can be empty, in which case, we keep all atoms (and not none). Example: {"T": {"min": 300, "max": 350}}. The format of this dictionary must match the official format of an ROI (see the apply_ROI function for more details).
    ROD : dict
        Region of Disinterest: we exclude all atoms that are in this region. WARNING: this function acts axis by axis for now, unfortunately.
    find_arrival_times : bool
        If True, we fit the arrival time of the BEC. If False, we don't.
    n_max_cycles : int
        Maximum number of cycles we want to select. Not often useful.
    histogramm_width : float
        Width of the histogram bins for fitting according to T and .times, in ms.
    width_saturation : float
        Saturation width during which there is no signal due to TDC saturation. The histogram points between tmax, the time at which the signal is maximal, and tmax + dt are removed and not taken into account in the fit.
    supplementary_rois : list of ROIs
        Additional ROIs in which we want to count the atoms (to apply filters during sequence analysis, for example). Data from these ROIs will be added to the dataframe BEC arrival_times.
    metadata : list of string
        List of metadata to load with the sequence parameters.

    """
    ### STEP 1 : gather data and save it
    selected_files = select_atoms_in_folder(folder)
    ROI_for_fit = check_roi_for_fit(ROI_for_fit)
    N_files = min([len(selected_files), n_max_cycles])
    # on va remplir les listes suivantes avec l'ensemble des atomes.
    Xc, Yc, Tc, Cyclec = [], [], [], []
    selected_files = selected_files[0:N_files]
    print(f"Starting to gather atoms from {folder}")

    df_atoms = pd.DataFrame()
    df_parameters = pd.DataFrame()
    df_arrival_times = pd.DataFrame()
    ##  take off warnings
    warnings.simplefilter("error", OptimizeWarning)
    for i, filename in tqdm(enumerate(selected_files)):
        cycle_prefix = filename.replace(".atoms", "")
        X, Y, T = load_XYTTraw(filename)
        seq, cycle = return_cycle_from_path(filename, show_error=False)
        raw_data = pd.DataFrame({"X": X, "Y": Y, "T": T})
        raw_data["Cycle"] = cycle
        atoms_in_ROI = apply_ROI(raw_data, ROI)
        atoms_in_ROI = apply_ROD(atoms_in_ROI, ROD)
        df_atoms = pd.concat([df_atoms, atoms_in_ROI])

        ## On charge maintenant les parametres de séquence
        seq_par = load_metadata(cycle_prefix, "parameters")

        for meta in metadata:
            seq_par.update(load_metadata(cycle_prefix, meta, show_error=False))
        seq_par["Sequence"] = seq
        seq_par["Cycle"] = cycle
        seq_par["Cycle time"] = get_creation_time(filename)
        paris_timezone = pytz.timezone("Europe/Paris")
        seq_par["Date"] = datetime.fromtimestamp(
            get_creation_time(filename), tz=paris_timezone
        )
        # add ROI nuumber of atoms and so on
        seq_par["Number of Atoms in ROI"] = len(atoms_in_ROI)
        seq_par["dT of atoms in ROI"] = np.std(atoms_in_ROI["T"])
        for i, roi in enumerate(supplementary_rois):
            at = apply_ROI(raw_data, roi)
            seq_par[f"Number of Atoms in supplementary ROI{i}"] = len(at)
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
            df_arrival_times = pd.concat(
                [df_arrival_times, pd.DataFrame(arrival, index=[i])]
            )
    seq_id_str, seq_id = get_sequence_id_from_seq_dir(folder)
    df_parameters["Sequence ID"] = seq_id
    df_parameters["Sequence ID str"] = seq_id_str
    if not find_arrival_times and "center (ms)" in df_parameters.keys():
        df_parameters["BEC Arrival Time"] = df_parameters["center (ms)"]
    filename_dataset = os.path.join(folder, "dataset.pkl")
    df_atoms.to_pickle(filename_dataset)
    filename_parameters = os.path.join(folder, "parameters.pkl")
    df_parameters.to_pickle(filename_parameters)
    if find_arrival_times:
        filename_arrival_time = os.path.join(folder, "arrival_times.pkl")
        df_arrival_times = df_arrival_times.merge(df_parameters, on="Cycle")
        df_arrival_times.to_pickle(filename_arrival_time)
        return [filename_dataset, filename_parameters, filename_arrival_time]

    return [filename_dataset, filename_parameters, None]


def export_metadatas_to_pickle(folder, metadata_list=["json"]) -> pd.DataFrame:
    """Export all metadatas of a folder to a pickle file. Return a dataframe with
    all metadatas in the metadata_list.

    Parameters
    ----------
    folder : str
        path of the folder from wichi you want to load metadatas.
    metadata_list : list, optional
        list of string with metadatas you want to load. By default ["json"]

    Returns
    -------
    pdDataFrame
        dataframe containg all metadatas from the metadata_list
    """
    all_cycles = glob.glob(os.path.join(folder, "*.sequence_parameters"))
    df = pd.DataFrame()
    for i, cycle_path in tqdm(enumerate(all_cycles)):
        cycle_prefix = cycle_path.replace(".sequence_parameters", "")
        seq, cycle = return_cycle_from_path(cycle_prefix)
        seq_par = {}
        for metadata in metadata_list:
            seq_par.update(load_metadata(cycle_prefix, metadata))
        seq_par["Sequence"] = seq
        seq_par["Cycle"] = cycle
        seq_par["Cycle time"] = get_creation_time(cycle_path)
        paris_timezone = pytz.timezone("Europe/Paris")
        seq_par["Date"] = datetime.fromtimestamp(
            get_creation_time(cycle_path), tz=paris_timezone
        )
        new_df = pd.DataFrame(seq_par, index=[i])
        df = pd.concat([df, new_df])
    seq_id_str, seq_id = get_sequence_id_from_seq_dir(folder)
    df["Sequence ID"] = seq_id
    df["Sequence ID str"] = seq_id_str
    df = df.reset_index(drop=True)
    df.to_pickle(os.path.join(folder, "metadata.pkl"))
    return df


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
                key = element["name"]
                if element["unit"] != "":
                    key += "  ({})".format(element["unit"])
                value = element["value"]
                seq_par[key] = value

    return seq_par


def gather_saved_sequence_parameters(folder):
    """This function retrieves all the parameters saved for each cycle of the sequence and puts them into a dataframe with the cycle number

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
    filename_parameters = os.path.join(folder, "parameters.pkl")
    dataframe.to_pickle(filename_parameters)
    return dataframe


def load_metadata(cycle_prefix, metadata, show_error=True):
    """load the metadata dictionary associated to the cycle prefix.

    Parameters
    ----------
    cycle_prefix : string
        prefix of the cycle. For example /mnt/...../002_005
    metadata : string
        type of metadata to load.
        For example 'fit' for HAL fit or 'pico' for picoscope files.

    Returns
    -------
    dict
        dictionary with associated metadatas
    """
    if metadata.lower() in ".json parameters params":
        return load_hal_type_metadata(cycle_prefix + ".json", show_error=show_error)
    if metadata.lower() in ".picoscope_treated":
        return load_hal_type_metadata(
            cycle_prefix + ".picoscope_treated", show_error=show_error
        )
    if metadata.lower() in ".picoscope_treated2000 arrival time bec arrival time":
        return load_hal_type_metadata(
            cycle_prefix + ".picoscope2000_treated", show_error=show_error
        )
    if metadata.lower() in "HAL fits camera .HAL_fits":
        directory, file_name = os.path.split(cycle_prefix)
        file_name = file_name + ".json"
        file_name = os.path.join(
            directory, ".HAL_fits", file_name, show_error=show_error
        )
        return load_hal_type_metadata(file_name, show_error=show_error)
    if metadata.lower() in "all parameters every parameters":
        return load_dictionary_metadata(
            cycle_prefix + ".sequence_parameters", show_error=show_error
        )
    if metadata.lower() in ".mcp stats .mcp_stats":
        directory, file_name = os.path.split(cycle_prefix)
        file_name = file_name + ".json"
        file_name = os.path.join(directory, ".MCPstats", file_name)
        return load_hal_type_metadata(file_name, show_error=show_error)
    return {}


import os


def load_hal_type_metadata(file, show_error=True) -> dict:
    """Load HAL type of metadat i.e. a list of dictionary with
    entries name, value, unit, error etc...

    Args:
        file (str or pathlib.Path): path to the metadata file
        show_error (boolean): if we show the log

    Returns:
        dict: dictionary containing all HAL metadata of the file.
    """
    try:
        dico = {}

        f = open(file)
        data = json.load(f)

        for element in data:
            key = element["name"]
            if element["unit"] != "":
                key += " ({})".format(element["unit"])
            dico[key] = element["value"]
        return dico
    except Exception as e:
        msg = f"{__file__}"
        msg += " \n     from load_hal_type_metadata \n "
        msg += f"Loading HAL type file {file} failed. Are you sure "
        msg += f"the file you want to load is a HAL file ? Error is {e}."
        if show_error:
            log.error(msg)
        return {}


def load_dictionary_metadata(file, key_separator=" | ", show_error=True) -> dict:
    """load a dictionary from file, flatten it with separator and returns

    Parameters
    ----------
    file : string or pathlib path
        path to the file you want to load
    key_separator : str, optional
        _description_, by default " | "

    Returns
    -------
    dict
        flatten dictionary from file
    """
    try:
        f = open(file)
        data = json.load(f)
        reducer = reducers.make_reducer(delimiter=key_separator)
        data = flatten(data, reducer=reducer)
        return data
    except Exception as e:
        msg = f"{__file__}"
        msg += " \n     from load_dictionary_metadata \n "
        msg += f"Loading dictionary from {file} failed. Are you sure "
        msg += f"the file you want to load is a dictionnary-like file ? Error is {e}."
        if show_error:
            log.error(msg)
    return {}


def get_creation_time(file_path):
    # Obtient les métadonnées du fichier
    file_stat = os.stat(file_path)

    # Choisis la bonne index pour l'heure de création en fonction du système d'exploitation
    if platform.system() == "Windows":
        creation_time = file_stat.st_ctime
    else:
        # Sur les systèmes Unix, st_birthtime est l'heure de création (si disponible)
        # Sinon, utilise st_mtime (dernière modification)
        creation_time = getattr(file_stat, "st_birthtime", file_stat.st_mtime)

    return creation_time


if __name__ == "__main__":
    # folder = Path("/mnt/manip_E/2022/11/18/083")
    # export_data_set_to_pickle(
    #     folder,
    #     ROI={"T": {"min": 306.2, "max": 309.7}},
    #     find_arrival_times=True,
    #     n_max_cycles=3,
    # )
    if False:
        cycle_prefix = "/mnt/manip_E/2023/10/09/003/003_005"
        dico = load_metadata(cycle_prefix, "mcp")
        print(dico)

    if True:
        root = "/mnt/manip_E/2023/10"
        all_seq = [
            "05/007",
            "05/011",
            "05/013",
            "05/016",
            "05/017",
            "06/009",
            "06/013",
            "07/017",
            "07/019",
            "07/021",
            "07/025",
            "08/001",
            "08/004",
        ]
        all_seq = ["08/011", "09/001"]
        for seq in all_seq:
            folder = os.path.join(root, seq)
            export_metadatas_to_pickle(folder, ["pico", "all parameters"])
    if False:  # check du fit du BEC
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
        folder = "/mnt/manip_E/2023/09/22/060"
        export_data_set_to_pickle(
            folder,
            ROI={"T": [312, 322]},
            find_arrival_times=False,
            metadata=["picoscope"],
        )
        data = pd.read_pickle("/mnt/manip_E/2023/09/22/060/parameters.pkl")
        print(data.columns)
        import matplotlib.pyplot as plt

        plt.plot(data["Cycle"], data["Bragg fit error"])
        plt.show()
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
