import pandas as pd
import logging
from heliumtools.misc.logger import getLogger
import numpy as np

log = getLogger(__name__)
# set the desired level of warning : DEBUG / INFO / WARNING
log.setLevel(logging.INFO)


def data_filter(data, bec_arrival_times, filters):
    """filters the dataframe bec_arrival_time based on the filters. It then filter out the "Cycles" of data which are NOT anymore in bec_arrival_time.
    Exemple: data_filter(data, bec_arrival_times, {"BEC Arrival time":[307, 308]}) will select all cycles in bec_arrival_times which are within [307,308] ms. We then select only the Cycles of data which have the cycles present in bec_arrival_time.

    Parameters
    ----------
    data : pandas.DataFrame
        A DataFrame containing a columns "Cycle" to compare with bec_arrival_time (often an 'atom' dataframe)
    bec_arrival_times : pandas.DataFrame
        _description_
    filters : dic
        ROI type dictionary e.g. {"key":[min, max]}, or {"key":{"center":mean, "size":size}} and so on....

    Returns
    -------
    (pandas.DataFrame, pandas.DataFrame)
        data and bec_arrival_times dataframe that have been filtered out.
    """
    data = data.reset_index(drop=True)
    bec_arrival_times = bec_arrival_times.reset_index(drop=True)
    selec_bec_arrival_times = apply_ROI(bec_arrival_times, filters)

    selected_data = data[data["Cycle"].isin(selec_bec_arrival_times["Cycle"])]
    return selected_data, selec_bec_arrival_times


def apply_ROD(df, ROD):
    """This function returns athe atoms dataframe such that all elements are OUTSIDE the ROD dictionary range. If the ROD is an empty dictionary, this function returns the initial dataframe.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame with all atoms.
    ROI : dic
        dictionary for which every entry (for exemple 'T') matches a column of the df dataframe, e.g. {"key":[min, max]}, or {"key":{"center":mean, "size":size}} and so on....

    Returns
    -------
    pandas.DataFrame
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
            log.warning(
                f"[heliumtools.tools.apply_ROD] The key {key} of the ROD is not in the other dataframe."
            )
    return df


def apply_ROI(atoms, ROI):
    """
    This function returns the atoms dataframe such that all elements are within the ROI dictionary range. If the ROI is an empty dictionary, this function returns the initial dataframe.

    Parameters
    ----------
    atoms : pandas dataframe
        DataFrame with all atoms.
    ROI : dic
        dictionary for which every entry (for exemple 'T') matches a column of the atoms dataframe. e.g. {"key":[min, max]}, or {"key":{"center":mean, "size":size}} and so on....

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
                log.warning(
                    f"[heliumtools.tools.apply_ROI]: the key ({key}) of your ROI is not in the dataframe on which you want to apply ROI."
                )

    return atoms


def get_roi_min_max(roi, key):
    """Base function that returns the maximum and minimum value of a roi type dictionary for a given entry key. This function is called by the apply_ROI and ROD methods.
    Current implemented ROI format:
        {"key":{"range":[-3,3]}}
        {"key":{"minimum":-3, "maximum":3}}
        {"key":{"min":-3, "max":3}}
        {"key":{"center":0, "size":6}}
        {"key":{"position":0, "size":6}}
        {"key":[-3,3]}
        {"key":(-3,3)}

    Parameters
    ----------
    roi : dictionary
        ROI type dictionary
    key : str
        key of the dictionary from which we we want the minimum and maximum value

    Returns
    -------
    tuples
        maximum and minimum value of the key of the roi.
    """
    if key not in roi:
        log.warning(
            f"[heliumtools.tools.get_roi_min_max] the axis {key} is not in the {ROI}."
        )
        return (-np.inf, np.inf)
    value = roi[key]
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
    elif "center" in value and "size" in value:
        minimum = value["center"] - 0.5 * value["size"]
        maximum = value["center"] + 0.5 * value["size"]
    elif type(value) == list or type(value) == tuple or type(value) == np.ndarray:
        minimum = min(value)
        maximum = max(value)
    else:
        log.warning(
            "[heliumtools.tools.get_roi_min_max] The ROI format of {} at the key {key} was not recognized. Please read the get_roi_min_max documentation. We expect a dictionary with all values being either a dictionary or a list. ".format(
                roi, key
            )
        )
        return (-np.inf, np.inf)
    return (minimum, maximum)


def get_roi_size(roi, axis) -> float:
    """Return the size of a roi along a given axis.

    Parameters
    ----------
    roi : dictionary
        ROI type dictionary
    key : str
        key of the dictionary from which we we want the minimum and maximum value

    Returns
    -------
    float
        Returns the size of a ROI like dictionary along a given axis.
    """
    (minimum, maximum) = get_roi_min_max(roi, axis)
    return maximum - minimum


def get_roi_center(roi, axis) -> float:
    """Return the center of a roi along a given axis.

    Parameters
    ----------
    roi : dictionary
        ROI type dictionary
    key : str
        key of the dictionary from which we we want the minimum and maximum value

    Returns
    -------
    float
        Returns the center of a ROI like dictionary along a given axis.
    """
    (minimum, maximum) = get_roi_min_max(roi, axis)
    return (maximum + minimum) / 2


def bootstrap_dataframe(original_data, key="Cycle"):
    """Bootstrap a dataframe if it is regular (i.e. the numbe rof row per cycle is alays the same).
    See 03/01/2024 for motivation.

    Parameters
    ----------
    original_data : dataframe
        _description_
    key : str, optional
        clé du dataframe sur lequel on veut faire le bootstrap, by default "Cycle"

    Returns
    -------
    pandas dataframe
        the new data dataframe bootstrapped
    """
    initial_keys = original_data[key].unique()
    n_cycles = len(initial_keys)
    original_data["Initial " + key] = original_data[key]
    ordata = original_data.set_index([key, original_data.groupby(key).cumcount()])
    # ordata.index.names = [key, "My_tmp_index"]
    original_data_array = ordata.values.reshape((n_cycles, -1, len(ordata.columns)))
    _, n_tmp_index, _ = original_data_array.shape
    new_indices = np.random.randint(0, n_cycles, n_cycles)
    NEW = np.zeros_like(original_data_array)
    NEW[:] = original_data_array[new_indices]
    data = pd.DataFrame(
        data=NEW.reshape((-1, len(ordata.columns))), columns=list(ordata.columns)
    )
    data["Cycle"] = np.repeat(initial_keys, n_tmp_index)
    # data["N"] = np.tile(np.arange(n_tmp_index), n_cycles)
    return data
