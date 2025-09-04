import numpy as np

# define square pulse
def SquarePulse(t,t1,t2):
    """ Square pulse 
    Parameters
    ---------------------
    t : numpy array
        time
    t1 : float
        begining of pulse
    t2 : float
        end of pulse
    Return
    -----------------------
    Pulse value
    """
    pulse = np.zeros(len(t),dtype=float)
    mask = (t >= t1)*(t <= t2)
    pulse[mask] = 1.0
    return pulse

# define sinc pulse
def SincPulse(t,t1,t2,OmegaM,splitter = True):
    """ Sinc pulse 
    Parameters
    ---------------------
    t : numpy array
        time
    t1 : float
        begining of pulse
    t2 : float
        end of pulse
    OmegaM : float
        2*pi*(Rabi frequency). It should have same units as t 
    splitter : bool
        True for sinc in splitter configuration. False otherwise.
    Return
    -----------------------
    Pulse value
    """
    # compute sinc OmegaS
    if splitter:
        OmegaS = 2*OmegaM
    else:
        OmegaS = OmegaM
    
    # pulse duration
    tau = t2 - t1
    # compute pulse
    return np.sinc(OmegaS*(t-tau/2-t1)/np.pi)*SquarePulse(t,t1,t2)

# define reburp pulse
def ReburpPulse(t,t1,OmegaM):
    """ Square pulse 
    Parameters
    ---------------------
    t : numpy array
        time
    t1 : float
        begining of pulse
    OmegaM : float
        2*pi*(Rabi frequency). It should have same units as t 
    Return
    -----------------------
    Pulse value
    """
    # Fourier coefficients An for the reburp pulse
    reburp_coefficients = [0.48, -1.03, 1.09, -1.59, 0.86, -0.44, 0.27, -0.17, 0.10, -0.08, 0.04, -0.04, 0.01, -0.02, 0.00, -0.02]
    orders = np.arange(len(reburp_coefficients))
    
    # compute reburp Omega_S parameter
    OmegaS = 2*reburp_coefficients[0]*OmegaM
    # compute pulse duration
    tau = 2*np.pi/OmegaS

    # pulse
    pulse = np.zeros(len(t))
    # define mask for pulse on
    mask = (t >= t1)*(t <= (tau + t1))
    #  calculate series
    arg = np.tensordot(orders*OmegaS,t[mask]-t1,axes = 0)
    cossines = np.cos(arg)
    pulse[mask] = np.tensordot(reburp_coefficients,cossines,axes = 1) 

    return pulse
    


