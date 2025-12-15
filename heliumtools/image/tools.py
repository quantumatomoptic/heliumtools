from PIL import Image
import numpy as np
import scipy.optimize as opt
import os , glob , json
from flatten_dict import flatten, reducers
from tqdm import tqdm

""" Define functions to load and manipulate images """

# Transforms camera image into a 2D numpy array
def ImageToArray(image):
    """ Function that given the path to a png image with grey scale format, it transforms it into a 2D numpy array
    Parameters
    ----------------------------------------
    image : str
        path ot image
    Return
    ------------------------------------------
    2D numpy array format of image.
    """
    # load image
    img = Image.open(image)
    # return 2D numpy array
    return np.asarray(img)

# Fits a 2D numpy array
def FitImage(func,x,y,image,guess):
    """ Function that fits a 2D numpy array that is the image.
    Parameters
    ----------------------------------------------------------
    image : 2D numpy array
        2D numpy array of image
    x : 1D numpy array
        x-axis values of image
    y : 1D numpy array
        y-axis values of image
    func : python function
        function to fit image. Function should have as arguments func(XY,*args) where XY should of the format [x,y] where x and y are the x and y axis values.
        args are the function fit parameters
    guess : numpy array
        array with guess of fit parameters
    Return
    ---------------------------------------------------------------
    popt : optimized parameters
    pcov : covariance matrix
    """

    # create meshgrid
    X , Y = np.meshgrid(x,y)
    # get shape of grid
    size = X.shape
    # resize X and Y to flatten them
    x_1d = X.reshape((1, np.prod(size)))
    y_1d = Y.reshape((1,np.prod(size)))
    # stack axis
    xy_stack = np.vstack((x_1d, y_1d))
    # flatten image array for fit
    image = image.flatten()

    # try to fit
    try:
        popt , pcov = opt.curve_fit(func,xy_stack,image,p0=guess)
    except:
        print("Fit failed!")
        popt = guess
        pcov = []

    return popt , pcov
   
# Applies a Roi to a 2D numpy array
def ApplyRoi(image,x,y,Roi):
    """ Applies a Roi to the image.
    Parameters
    -----------------------------------------------
    image : 2D numpy array
        2D numpy array of image
    x : 1D numpy array
        x-axis coordinates
    y : 1D numpy array
        y-axis coordinates
    Roi : dictionary of the type {"x" : [min,max]  , "y" : [min,max]} with min and max values for each axis
    Return 
    --------------------------------------------------
    Cropped image 2D numpy array
    New x and y axis values
    """

    # compute x and y mask
    Xmask = (x >= Roi["x"][0])*(x <= Roi["x"][1])
    Ymask = (y >= Roi["y"][0])*(y <= Roi["y"][1])

    # apply masks to axis
    x = x[Xmask]
    y = y[Ymask]

    # create mesh of masks
    Xmask, Ymask = np.meshgrid(Xmask,Ymask)
    # multiply X and Y masks and flatten it
    mask = (Xmask*Ymask).flatten()
    # apply mask to flattened image
    image  = image.flatten()[mask]

    # return reshaped image
    return image.reshape((y.shape[0],x.shape[0])) , x , y 

# gets all filepaths for images
def gatherImages(folder,sequences):
    """ Gets all images filepaths.
    Parameters
    --------------------------------------------------------
    folder : str
        path to folder where sequences with images are. E.g "/mnt/manip_E/2025/05"
    sequences : list of str
        list with all sequences to gather
    Returns 
    ----------------------------------------------------------
    List with all filepaths
    """
    
    extension = "*.png"
    
    # Iterate over each sequence and collect filepaths
    files = []
    for seq in sequences:
        filepaths = glob.glob(os.path.join(folder, seq, extension))
        files.extend(filepaths)
    
    return files

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


""" Define usefull fit functions for imaging """

def Gaussian2D(XY,x0,y0,A0,sigmaX,sigmaY,CovXY,offset):
    """ 2D gaussian function.
    Parameters
    --------------
    XY : numpy array
        axis values. XY should of the format [x,y] where x and y are the x and y axis values.
    x0 , y0 : floats
        coordinates of the center of the distribuition
    A0 : float
        center of the distribuition
    sigmaX , sigmaY , CobXY : floats
        Covariance matrix elements
    offset : float
        offset of the distribuition
    Returns
    --------------
        Skewed 2D gaussian function
    """
    x = XY[0]-x0
    y = XY[1]-y0
    # calculate gaussian function with covariance
    gauss = np.exp(-np.power(x/sigmaX,2)/2)*np.exp(-np.power(y/sigmaY,2)/2)*np.exp(-CovXY*x*y)
    return A0*gauss + np.abs(offset)
