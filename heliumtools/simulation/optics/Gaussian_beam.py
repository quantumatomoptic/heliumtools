import numpy as np

""" This is a library with useful optics functions for gaussian beam calculations """

""" define some useful classes """

# define class to compute waist of the beam before and after lens
class WaistClass():
    """ Class to compute waist of the beam before and after a lens transformation.
    Parameters
    ------------------------------------
    wavelength : float
        laser's wavelength
    """

    def __init__(self,wavelength):
        self.wavelength = wavelength

    def LensSolver(self,center,f,w0,LensPosition):
       """ Computes the transformation of a Gaussian beam by a thin lens
        Parameters
        --------------------------------------
        center : float
            position of the waist
        f : float
            focal length of the lens
        w0 : float
            original waist 1/e^2 radius 
        LensPosition : float
            position of the lens
        Return
        --------------------------------------
        newBeam : dict
            dictionary with new Gaussian beam parameters {"M" , "w0" , "waist position" , "z0"}
        """
       # compute lens transformation. LensSolver() assumes lens at origin
       newBeam = LensSolver(LensPosition - center,f,w0,self.wavelength)
       # transform distance to coordinate
       newBeam["waist position"] = newBeam["waist position"] + LensPosition
       return newBeam

    def waist(self,z,w1,P1,f,LensPosition):
        """ computes waist of the beam at position z 
        Parameters
        ------------------------------------------------
        z : numpy array
            position relative to the origin at which we compute waist
        w1 : float
            original 1/e^2 radius waist before lens transformation
        P1 : float
            position of the beam waist relative to the origin before lens transformation
        f : float
            focal length of the lens
        LensPosition : float
            position of the lens relative to the origin
        Return 
        -----------------------------------------------
            waist at position z
        """

        # initialize waist array
        w = np.zeros(len(z),dtype = float)

        """ if we are before the lens """
        mask = z < LensPosition
        w[mask] = WaistFunc(z[mask],w1,P1,self.wavelength)

        """ if we are after the lens """
        # compute new waist and its position
        newBeam = self.LensSolver(P1,f,w1,LensPosition)
        # compute waist transformation
        w2 = newBeam["w0"] # new beam waist
        P2 = newBeam["waist position"] # new beam waist position

        # compute waist
        mask = z > LensPosition
        w[mask] = WaistFunc(z[mask],w2,P2,self.wavelength)

        return w


# define more general class than WaistClass() to compute Guassian beam under general transformations, e.g ABCD matrices
class BeamClass():
    """ More general class than WaistClass() to compute Guassian beam under general transformations, e.g ABCD matrices.
    Parameters
    ------------------------------------
    wavelength : float
        laser's wavelength
    """
    def __init__(self,wavelength):
        self.wavelength = wavelength

    # initiate q at waist position
    def get_q0(self,w0):
        """ initiate beam q(z) function at waist position """
        z0 = Rayleigh_length(w0,self.wavelength)
        return 1j*z0
    
    # compute q(z) given initial beam waist w1 and its center P1
    def getQ(self,z,w1,P1,LensPosition,Mlens):
        """ Computes q(z) given initial beam waist w1 and its center P1
        Parameters
        -----------------------------------------------------------
        z : numpy array
            position relative to the origin at which we compute waist
        w1 : float
            original 1/e^2 radius waist before lens transformation
        P1 : float
            position of the beam waist relative to the origin before lens transformation
        f : float
            focal length of the lens
        LensPosition : float
            position of the lens relative to the origin
        Mlens : 2D matrix
            ABCD matrix of a general optics element, e.g a lens
        Return 
        -----------------------------------------------
            q(z) parameter at z
        """
        # compute q(z) at beam waist
        q1 = self.get_q0(w1)

        # propgate up to lens
        q2 = prodLaw(q1,air(LensPosition-P1))
        # apply lens ABCD matrix
        q2 = prodLaw(q2,Mlens)

        # initialize Q function
        Q = np.zeros(len(z),dtype = complex)

        # for each value of z
        for i in range(0,len(z)):
            """ if we are before the lens """
            if z[i] <= LensPosition:
                Q[i] = prodLaw(q1,air(z[i]-P1))
            
            """ if we are after the lens """
            if z[i] > LensPosition:                
                Q[i] = prodLaw(q2,air(z[i]-LensPosition))

        return Q
    
    # compute waist of beam at z from q(z) parameter
    def waist(self,z,w1,P1,LensPosition,Mlens):
        """ Computes waist of beam at z from q(z) parameter
        Parameters
        -----------------------------------------------------------
        z : numpy array
            position relative to the origin at which we compute waist
        w1 : float
            original 1/e^2 radius waist before lens transformation
        P1 : float
            position of the beam waist relative to the origin before lens transformation
        f : float
            focal length of the lens
        LensPosition : float
            position of the lens relative to the origin
        Mlens : 2D matrix
            ABCD matrix of a general optics element, e.g a lens
        Return 
        -----------------------------------------------
            waist at z
        """

        # compute beam q(z) function
        Q = self.getQ(z,w1,P1,LensPosition,Mlens)
        # compute waist
        return getWaist(Q,self.wavelength)
    
    def curvature(self,z,w1,P1,LensPosition,Mlens):
        """ Computes curvature of beam at z from q(z) parameter
        Parameters
        -----------------------------------------------------------
        z : numpy array
            position relative to the origin at which we compute waist
        w1 : float
            original 1/e^2 radius waist before lens transformation
        P1 : float
            position of the beam waist relative to the origin before lens transformation
        f : float
            focal length of the lens
        LensPosition : float
            position of the lens relative to the origin
        Mlens : 2D matrix
            ABCD matrix of a general optics element, e.g a lens
        Return 
        -----------------------------------------------
            waist at z
        """
        # compute beam q(z) function
        Q = self.getQ(z,w1,P1,LensPosition,Mlens)
        # compute waist
        return getCurvature(Q,self.wavelength)


""" define some useful functions to compute waist and lens transformation for a Gaussian beam """

# computes the transformation of the beam by a thin lens
def LensSolver(z,f,w0,wavelength):
    """ Computes the transformation of a Gaussian beam by a thin lens
    Parameters
    --------------------------------------
    z : float or numpy array
        the distance of the lens relative to the original beam center. z should be positive if lens is after the waist and negative otherwise.
    f : float
        focal length of the lens
    w0 : float
        original waist 1/e^2 radius 
    wavelength : float
        laser's wavelength
    Return
    --------------------------------------
    newBeam : dict
        dictionary with new Gaussian beam parameters {"M" , "w0" , "waist position" , "z0"}
    """
    

    # compute Rayleigh length 
    z0 = Rayleigh_length(w0,wavelength)

    # dictionary to hold new beam parameters
    newBeam = dict()

    # compute magnification
    newBeam["M"]  = np.abs(f)/np.sqrt(np.power(z-f,2)+z0**2)

    # compute new waist
    newBeam["w0"] =  newBeam["M"]*w0

    # compute waist position. Distance is from the lens.
    # if the results is positive, then the new waist position is after the lens, otherwise it is negative
    newBeam["waist position"] =  f + np.power(newBeam["M"],2)*(z-f)

    # compute new Rayleigh distance
    newBeam["z0"] =  np.power(newBeam["M"],2)*z0

    return newBeam

# computes beam waist at position z
def WaistFunc(z,w0,center,wavelength):
    """ Computes beam waist at position z
    Parameters
    ------------------------------------------------------
    z : float or numpy array
        position at which we compute the waist
    center : float
        position of the waist
    w0 : float 
        waist 1/e^2 radius 
    wavelength : float
        laser's wavelength
    Return 
    ------------------------------------------------------
        waist at position z
    """

    # compute Rayleigh length 
    z0 = Rayleigh_length(w0,wavelength)
    # compute function
    return w0*np.sqrt(1 + np.power((z-center)/z0,2))

# computes Rayleigh length 
def Rayleigh_length(w0,wavelength):
    return np.pi*np.power(w0,2)/wavelength


""" define ABDC matrices for propagation of a beam """

# matrix for lens transfer
def lens(f):
    return np.array([[1,0],[-1/f,1]],dtype = complex)

# matrix for propagation of beam though air or any material
def air(d):
    return np.array([[1,d],[0,1]],dtype = complex)

# def refraction at surface with curvature R matrix
def refraction(R,n1,n2):
    return np.array([[1,0],[-(n2-n1)/(n2*R),n1/n2]],dtype = complex)

# def refraction at plane surface
def refraction_plane(n1,n2):
    return np.array([[1,0],[0,n1/n2]],dtype = complex)

def prodLaw(q,M):
    return (M[0,0]*q + M[0,1])/(M[1,0]*q + M[1,1])

""" define useful functions to relate q(z) to R(z) and W(z) """

# get W(z) from q(z)
def getWaist(q,wavelength):
    W = np.sqrt(-wavelength/(np.pi*np.imag(1/q)))
    return W

# get R(z) from q(z)
def getCurvature(q):
    R = 1/np.real(1/q)
    return R 

# get q(z) given position z and Rayleigh range z0
def getQ(z,z0):
    return z + 1j*z0
