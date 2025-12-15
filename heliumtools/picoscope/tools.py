import os , glob , json
from flatten_dict import flatten, reducers
from tqdm import tqdm
import pandas as pd


# gets all filepaths for files with picoscope raw data
def gatherRawFiles(folder,sequences,picoscope):
    """ Gets all filepaths for files with picoscope raw data.
    Parameters
    --------------------------------------------------------
    folder : str
        path to folder where sequences with data are. E.g "/mnt/manip_E/2025/05"
    sequences : list of str
        list with all sequences to gather
    picoscope : str
        from which picoscope you want data. Implemented options are "picoscope 3000" and "picoscope 2000" 
    Returns 
    ----------------------------------------------------------
    List with all filepaths
    """
    # if want the data from picoscope 3000
    if picoscope == "picoscope 3000":
        extension = "*.picoscope_raw_data"
    # if we want the data from picosope 2000
    elif picoscope == "picoscope 2000":
        extension = "*.picoscope2000_raw_data"
    # else picoscope is not implemented
    else:
        print("Picoscope you ask for is not implemented!")
        return []

    # Iterate over each sequence and collect filepaths
    files = []
    for seq in sequences:
        filepaths = glob.glob(os.path.join(folder, seq, extension))
        files.extend(filepaths)
    
    return files

# function that read picoscope raw file
def readPicoRaw(file):
    """ function that read picoscope raw file.
    Parameters
    ------------------------------------------ 
    file : str
        path of file with raw data
    Returns
    -----------------------------------------
    data from file
    """
    return pd.read_pickle(file)

# loads metadata
def loadSeqParameters(folder,sequences):
    """ Gets sequence parameters.
    Parameters
    --------------------------------------------------------
    folder : str
        path to folder where sequences with data are. E.g "/mnt/manip_E/2025/05/23"
    sequences : list of str
        list with all sequences to gather
    Return 
    ---------------------------------------------------------
    Pandas dataframe with all sequence parameters
    """

    # sequence parameters file extension
    extension = "*.sequence_parameters"

    # dataframe to hold sequence parameters
    metadata = pd.DataFrame()

    # for each sequence
    counter = 0
    for seq in sequences:
        filepaths = glob.glob(os.path.join(folder, seq, extension))
        # for each file path
        print("Getting sequence parameters from: "+os.path.join(folder, seq))
        for file in tqdm(filepaths):
            # load dictionary
            df = load_dictionary_metadata(file)
            # add cycle id to dictionary
            df["cycle"] = int(df["cycle prefix"].split("/")[-1].split("_")[-1])
            df["cycle id"] = df["cycle"] + counter
            df["sequence number"] = int(df["sequence number"])
            # transform to pandas dataframe
            df = pd.DataFrame(df , index = [0])
            # apped to metadata
            metadata = pd.concat((metadata,df))
        # update counter for next sequence
        counter = counter + len(filepaths) + 1
    
    return metadata.sort_values(by = "cycle id").reset_index().drop(columns= ["sequence parameter path","scan parameter path","sequence folder","index"])


def load_dictionary_metadata(file, key_separator = " | ", show_error = True) -> dict:
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
