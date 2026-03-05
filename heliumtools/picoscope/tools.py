import os , glob , json
from flatten_dict import flatten, reducers
from tqdm import tqdm
import pandas as pd
from heliumtools.misc.gather_data import load_dictionary_metadata


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
    # if we want piscoscope 2000 phase
    elif picoscope == "picoscope 2000 phase":
        extension = "*.picoscope2000phase_raw_data"
    # else picoscope is not implemented
    else:
        print("Picoscope you ask for is not implemented!")
        return []

    # Iterate over each sequence and collect filepaths
    files = []
    for seq in sequences:
        filepaths = sorted(glob.glob(os.path.join(folder, seq, extension)))
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
            # add cycle number to dictionary
            df["cycle"] = int(df["cycle prefix"].split("/")[-1].split("_")[-1])
            # add cycle id
            df["cycle id"] = getCycleId(file)
            # add sequence number
            df["sequence number"] = int(df["sequence number"])
            # transform to pandas dataframe
            df = pd.DataFrame(df , index = [0])
            # apped to metadata
            metadata = pd.concat((metadata,df))
        # update counter for next sequence
        counter = counter + len(filepaths) + 1
    
    return metadata.sort_values(by = "cycle id").reset_index().drop(columns= ["sequence parameter path","scan parameter path","sequence folder","index"])


# given image file path it retrieves the cycle id in the json file
def getCycleId(file):
    """ Given image file path it retrieves the cycle id in the json file.
    ---------------------------------------
    Parameters
        file : str
        path to file image
    ---------------------------------------
    Return 
        cycle id from json file
    """
    # replace extension to json
    file = file.split(".")[0]+".json"
    # open json file
    with open(file, 'r') as f:
        array = json.load(f)
    # choose cycle id element
    array = array[2]
    return array["value"]

# Gets sequence parameters for a particular file image
def getSequenceParameters(file):
    """ Gets sequence parameters for a particular file image.
    ---------------------------------------
    Parameters
        file : str
        path to file image
    ---------------------------------------
    Return 
        dictionary with metadata
    """
    # replace extension
    file = file.split(".")[0]+".sequence_parameters"
    # return metadata
    return load_dictionary_metadata(file)
