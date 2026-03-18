"""
Author   : Rui Dias
Created  : 05/03/2026

Comments : Contains useful functions to analyse lattice pairs data, such as:
    Fit of pairs velocity
    Compute center of reference frame
    Compute diffraction profile for pairs
    etc.
"""

from scipy.special import erf
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
from tqdm import tqdm

from heliumtools.tools import apply_ROI
from heliumtools.simulation.Bragg.Bragg import Bragg
import heliumtools.simulation.Bragg.pulses as pulses

# computes the expected Bragg diffraction profile and compares it with pairs density
def diffraction_profile(Vz,density,center,simulation):
    """ Computes the expected Bragg diffraction profile and compares it with pairs density
    --------------------------
    Parameters
    Vz : numpy array
        Vz values in mm/s at which density is computed
    density : numpy array
        density of the pairs along Vz
    center : dict()
        center velocity of the pairs in mm/s. Should have format {"1" : 70 , "2" : 120}
    simulation : dict()
        dictionary with all parameters necessary to compute the diffraction profile.
        Keys are: ['vB (mm/s)','Pulse delay (ms)','Rabi frequency (kHz)','Pulse duration (ms)','Bragg detuning (kHz)','Pulse type','Sinc splitter']
    --------------------------
    Return
        fig , ax : plot objects
    """

    # get Bragg velocity and change it to mm/ms
    BraggVelocity = simulation["vB (mm/s)"]*1e-3

    # instanciate Bragg class
    bragg = Bragg(vB = BraggVelocity)

    # compute detuning to diffract pairs to +1 order
    chaves = list(center.keys())
    for key in chaves:
        detuning = bragg.frequencyOrderDifference(1,center[key]*1e-3)/(2*np.pi) # in kHz
        print("To diffract pair "+str(key)+" to +1 order, detuning should be "+str(np.round(detuning,2))+" kHz")

    print("\n")

    # compute detuning to diffract pairs to -1 order
    chaves = list(center.keys())
    for key in chaves:
        detuning = bragg.frequencyOrderDifference(0,center[key]*1e-3)/(2*np.pi) # in kHz
        print("To diffract pair "+str(key)+" to -1 order, detuning should be "+str(np.round(detuning,2))+" kHz")

    print("\n")


    """ define/save Bragg pulse parameters """

    # define pulse initial time in ms
    bragg.par["Pulse beggining"] = simulation["Pulse delay (ms)"]
    # define rabi frequency in kHz
    bragg.par["Rabi frequency"] = simulation["Rabi frequency (kHz)"]
    # define pulse duration in ms
    bragg.par["Pulse duration"] = simulation["Pulse duration (ms)"]
    # define Bragg detuning
    bragg.par["Detuning"] = simulation["Bragg detuning (kHz)"]
    # define number of diffraction orders
    bragg.par["diffraction orders"] = 2

    """ define computation parameters """
    Nmax = bragg.par["diffraction orders"]
    # define time step for differential equation
    bragg.par["time step solver"] = 1/(10*Nmax*bragg.par["Bragg recoil frequency"]/(2*np.pi)) # ms
    # define time step for propagation of wavefunction outside pulse
    bragg.par["time step propagator"] = 0.1
    # update Bragg parameters
    bragg.update_parameters()
    # define beggining of computation
    t0 = 0.0
    # define end of computation
    tfinal = bragg.par["Pulse beggining"] + bragg.par["Pulse duration"] + 0.1

    """ initialize pulse """
    # create time array to compute pulse
    time = np.arange(t0,tfinal+bragg.par["time step solver"],bragg.par["time step solver"])
    # OmegaM parameter
    OmegaM = 2*np.pi*bragg.par["Rabi frequency"]
    # beggining of pulse
    t1 = bragg.par["Pulse beggining"]
    # end of pulse
    t2 = bragg.par["Pulse beggining"] + bragg.par["Pulse duration"]

    if simulation["Pulse type"] == "reburp":
        # compute reburp pulse
        pulse = pulses.ReburpPulse(time,t1,OmegaM)
    elif simulation["Pulse type"] == "square":
        # compute square pulse
        pulse = pulses.SquarePulse(time,t1,t2)
    elif simulation["Pulse type"] == "sinc":
        pulse = pulses.SincPulse(time,t1,t2,OmegaM,splitter = simulation["Sinc splitter"])
    else:
        print("Invalid pulse type. I will do a reburp instead.")
        # compute reburp pulse
        pulse = pulses.ReburpPulse(time,t1,OmegaM)

    """ compute pulse phase """
    # We assume a constant detuning and a frequency slope to compensate gravity
    phase = pulses.phaseRamp(time,bragg.par["Detuning"],bragg.par["Frequency slope to compensate gravity"])

    """ initialize pulse shape and phase """
    bragg.initializePulse(time,pulse,phase)

    """ initial wavefunction """
    # create initial wavefunction, which will be filled with ones
    psi0 = np.ones(len(Vz),dtype=complex)
    # create initial 2D-array wavefunction at order zero for computation
    psi0 = bragg.createInitalPsi(Vz*1e-3,psi0,0)

    """ compute psi velocity """
    t , psi = bragg.compute_psiVelocity(t0,tfinal,Vz*1e-3,psi0,100)
    # compute norm of last time value of psi
    psi = np.real(psi[-1]*np.conjugate(psi[-1]))

    """ plot """
    fig , ax = plt.subplots(ncols = 2 , figsize = (2*6,5))
    # plot -2 and -1 orders of diffraction
    ax[0].plot(Vz,psi[Nmax - 2],color = "black" , label = "n = -2")
    ax[0].plot(Vz,psi[Nmax - 1],color = "red", label = "n = -1")
    # plot +1 and +2 orders of diffraction
    ax[1].plot(Vz,psi[Nmax + 1],color = "darkgreen" , label = "n = 1")
    ax[1].plot(Vz,psi[Nmax + 2],color = "purple", label = "n = 2")
    # plot labels and density
    ax[0].set_ylabel("Diffraction efficiency",fontsize = 15)
    for i in [0,1]:
        ax[i].set_xlabel(r"$V_z$ (mm/s)",fontsize = 15)
        ax[i].legend(fontsize = 15)

        # plot horizontal lines
        ax[i].hlines(y = 1.0 , xmin = np.amin(Vz) , xmax = np.amax(Vz) , color = "grey" , ls = "dashed")
        ax[i].hlines(y = 0.5 , xmin = np.amin(Vz) , xmax = np.amax(Vz) , color = "grey" , ls = "dashed")

        # plot density
        ax2 = ax[i].twinx()
        color = 'tab:blue'
        if i == 1:
            ax2.set_ylabel('Pairs density', color=color,fontsize = 15)  # we already handled the x-label with ax1
        ax2.plot(Vz,density,color = color)
        ax2.tick_params(axis='y', labelcolor=color)

    return fig , ax

# computes density of pairs along Vz
def compute_density(atoms,Roi,binSize,step):
    """ Given a atoms dataframe with collumns Cycle, Vx, Vy and Vz, it computes the density of pairs along Vz.
    Notice that, the function allows for overlapping bins if step is smaller than binSize.
    ---------------------------------------
    Parameters:
        atoms : pandas dataframe
            dataframe with atoms collumns Cycle, Vx, Vy and Vz. E.g corr.atoms
        Roi : dictionary
            Roi for atoms.
        binSize : float
            Size of the bin of histogram
        step : float
            distance between bin centers
    ----------------------------------------
    Return:
        Vz : numpy array
            center of the bins
        density : numpy array
            histogram of the lattice pairs
    """

    # apply Roi to atoms
    atoms = apply_ROI(atoms,Roi)

    # get number of cycles
    Ncycles = len(atoms["Cycle"].drop_duplicates())

    # get min and max velocities
    min = np.amin(atoms["Vz"])
    max = np.amax(atoms["Vz"])

    # if step is equal to size of bin size, we use numpy histogram function
    if step == binSize:
        # create bins
        bins = np.arange(min,max+binSize,binSize)

        # build histogram
        density , Vz = np.histogram(atoms["Vz"], bins = bins)

        # normalize histogram and center bins
        density = density/Ncycles
        Vz = (Vz[1:]+Vz[0:-1])/2

    # if step is not equal to bin size
    else:
        # build bin centers
        Vz = np.arange(min,max+step,step)
        # recenter
        Vz = (Vz[1:]+Vz[0:-1])/2

        # initialize density array
        density = np.zeros(len(Vz),dtype = float)

        # count atoms in each bin
        for i in tqdm(range(0,len(Vz))):
            # define mask for bin
            mask = (atoms["Vz"] >= (Vz[i] - binSize/2))*(atoms["Vz"] < (Vz[i] + binSize/2))
            # count atoms
            density[i] = len(atoms["Vz"][mask])

        # normalize to number of cycles
        density = density/Ncycles


    return Vz , density

# fits pairs density and retrives the center of the pairs using a Skewed Gaussian function
def fit_pairs_density(Vz,density,FitRegion,guess,show):
    """ Given the density of the pairs along Vz, fits pairs density using a Skewed Gaussian function.
    The Skewed Gaussian takes 5 parameters: [center of the distribuition , amplitude , sigma , skeweness of gaussian , offset]. See SkewGaussian() for details
    --------------
    Parameters
    Vz : numpy array
        Vz values at which density is computed
    density : numpy array
        density of the pairs along Vz
    binSize : float
        size of bin to compute histogram of density
    FitRegion : dictionary
        Dictionay with fit regions for pair region 1 and pair region 2. It should have format:
            FitRegion = {"1" : {"min" : 45 , "max" : 85} , "2" : {"min" : 100 , "max" : 140} }
    guess : dictionary
        Fit guess parameters for skewed Gaussian for each region. It should have format:
            guess = {"1" : [70, 0.03,5, 0.5, 0] , "2" = [115,0.02,5,-0.5, 0] }
    show : bool
        True if you want to plot the fit result. False otherwise
    --------------
    Returns
        parameters : dictionary
            dictionary with fit parameters for each region
        center : dictionary
            center of pairs in each region
        CM : float
            shift in reference frame along Vz, such that pairs center of mass is at zero.
    """
    
    # get pairs labels
    pairs = list(FitRegion.keys())

    # if FitRegion and guess have different pairs labels
    if pairs != list(guess.keys()):
        print("FitRegion and guess dictionaries should have the same pairs labels.")
        return dict() , dict() , 0.0

    # create dataset to hold fit results
    dfFit = dict()
    for pair in  pairs:
        dfFit[pair] = []

    # create dataset to hold center of pairs values
    center = dict()
    for pair in pairs:
        center[pair] = 0.0

    """ fit pairs density in reach region, compute pairs center and center of mass """

    # initialize center of mass
    CM = 0.0

    # fit pair using skew gaussian
    for pair in pairs:
        # get density region inside Fit region 
        mask = (Vz >= FitRegion[pair]["min"])*(Vz <= FitRegion[pair]["max"])
        bins = Vz[mask]
        hist = density[mask]
        # try to fit
        try:
            popt, pcov = opt.curve_fit(SkewGaussian,xdata=bins,ydata=hist,p0=guess[pair])
        except:
            print("Fit of pair "+str(pair)+" failed !")
            popt = guess[pair]
            
        # add data to fit datasets
        dfFit[pair] = popt
        
        # find psition of skew gaussian where density is maximum
        center[pair] = findCenterSkewGaussian(popt)

        # compute center of mass
        CM = CM + center[pair]/2
        

    """ plot """
    if show:
        fig , ax = plt.subplots()
        # plot density
        ax.plot(Vz,density,color="tab:blue")
        # plot fit for each pair
        for pair in pairs:
            # get density region inside Fit region 
            mask = (Vz >= FitRegion[pair]["min"])*(Vz <= FitRegion[pair]["max"])
            bins = np.linspace(np.amin(Vz[mask]),np.amax(Vz[mask]),1000)
            # plot fit
            ax.plot(bins,SkewGaussian(bins,*dfFit[pair]),color="red")
            # plot center
            legenda = str(pair)+" center = "+str(np.round(center[pair],3))+" mm/s"
            ax.vlines(x = center[pair],ymin = 0,ymax = 1.1*np.amax(density),ls = "dashed",color = "black" , label = legenda)
        # plot center of mass
        ax.vlines(x = CM,ymin = 0,ymax = 1.1*np.amax(density),ls = "dashed",color = "green" , label = "CM")
        ax.set_xlabel(r"$V_z$ (mm/s)")
        ax.legend()
        plt.tight_layout()
        plt.show()


    return dfFit , center , CM

# define Skew gaussian to fil pairs
def SkewGaussian(x,x0,A0,sigma,alpha,offset):
    """ Skewed 2D gaussian function.
    --------------
    Parameters
    x : numpy array
        axis values
    x0 : floats
        coordinate of the center of the distribuition
    A0 : float
        amplitude of the distribuition
    sigma , : floats
        Covariance  element
    alpha : floats
        x  value
    offset : float
        offset of the distribuition
    --------------
    Returns
        Skewed 1D gaussian function
    """
    x = x-x0
    # calculate gaussian function with covariance
    gauss = np.exp(-np.power(x/sigma,2)/2)
    # calculate skewness 
    skewness = 1 + erf(alpha*x/np.sqrt(2))
    return A0*gauss*skewness + np.abs(offset)

# given a skew gaussian, finds where it is maximum
def findCenterSkewGaussian(args):
    """ Given a skew gaussian, finds where it is maximum. 
    -------------------------------------
    Parameters
        args : array
            array with Skewed Gaussian parameters. 
            Parameters should be [center,amplitude,sigma width,skewness,offset]
            See SkewGaussian() for details.
    -------------------------------------
    Return
        Position where Skewed gaussian is maximal
    """

    func = lambda x : -np.abs(SkewGaussian(x,*args))
    delta = args[2]/np.sqrt(1+args[2]**2)
    center = args[0]+args[1]*delta*np.sqrt(2/np.pi)
    return opt.fmin(func,x0=center,disp=False)[0]