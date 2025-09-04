import sys
import numpy as np
import scipy.integrate as integrate
from tqdm import tqdm, trange


import pulses as pulses

""" Units are mm/ms/kg """

class Bragg:
    
    """   Object initialization  """
    def __init__(self, **kwargs):

        """
        TODO:
            1) 
        """
        
        # dictionary to hold all relevant parameters in mm/ms/kg
        self.par = dict()

        """ Define constants in mm/ms/kg units """
        self.par["hbar"] = 1.054*10**(-31) # mm^2 kg/ms
        self.par["mHe"] = 6.65*10**(-27) # Helium's mass in kg
        self.par["gravity"] = 9.8067e-3 # mm/ms^2

        """ Define experimental lattice parameters in mm/ms/kg units """
        self.par["Bragg velocity"] = 50e-3 # Bragg velocity
        self.par["Rabi frequency"] = 1 # Rabi frequency
        self.par["Detuning"] = 10.2 # define detuning in kHz
        self.par["Phase slope"] = 0.0 # slope of frequency sweep in kHz/ms
        self.par["Pulse beggining"] = 0.0 # beggining of the pulse 
        self.par["Pulse duration"] = 1 # duration of the pulse
        self.par["Pulse type"] = "reburp" # type of pulse
        self.par["Splitter"] = True # In case you choose a sinc, you if want splitter configuration
        self.update_parameters() # update dictionary
        
        """ Define calculation parameters """
        self.par["velocity step"] = 1e-4 # define velocity step in mm/ms
        self.par["diffraction orders"] = 5 # define number of diffraction orders
        # define propagation step when pulse is not on
        self.par["time step propagator"] = 0.1
        # define time step to solve differential equation in ms when pulse is on
        if self.par["Rabi frequency"] > self.par["Bragg recoil frequency"]:
            self.par["time step solver"] = 1/(50*self.par["Rabi frequency"])
        else:
            self.par["time step solver"] = 1/(50*self.par["Bragg recoil frequency"])
        
        self.update_parameters() # update dictionary
        
    def update_parameters(self):
        """ Updates self.par which is a dictionary with all relevant experimental and calculation parameters """
        
        # update self.par
        self.par["mHe/hbar"] = self.par["mHe"]/self.par["hbar"]
        self.par["Bragg wavevector"] = self.par["mHe/hbar"]*self.par["Bragg velocity"]
        self.par["Bragg recoil frequency"] = np.power(self.par["Bragg wavevector"],2)/(2*self.par["mHe/hbar"])
        self.par["slope to compensate gravity"] = -self.par["Bragg wavevector"]*self.par["gravity"]
        self.par["Phase slope"] = self.par["slope to compensate gravity"]

    def frequencyOrderDifference(self,n,v):
        """ Frequency order difference (2 pi factor included) \delta_n = (\epsilon_{n+1} - \epsilon_{n})/hbar
        Parameters
        ------------------------
        n : float
            diffraction order to consider. E.g if n = 1, than FrequencyOrderDifference() gives the frequency difference between
            order n = 1 and order n = 0.
        v : float
            velocity of the atom in mm/ms units
        """
        delta0 = v*self.par["Bragg wavevector"] - self.par["Bragg recoil frequency"]
        return delta0 + 2*n*self.par["Bragg recoil frequency"]
   
    def braggPhase(self,t,detuning,slope):
        """ Time dependent phase of the beams. In the experiment we can only do a linear ramp, so the phase can only be quadratic.
        Parameters
        -----------------------------------
        t : numpy array or float
            time in ms
        detuning : float
            inital detuning between beams in kHz
        slope : float
            slope of the ramp detuning in kHz^2
        Return
        ----------------------------------------
        Phase at time t
        """
        return detuning*t + slope*np.power(t,2)/2

    def gpacket(self,t,t1,t2,OmegaM,type,splitter):
        """ Light pulse function.
        Parameters
        ---------------------
        t : float
            time in ms
        t1 : float
            begining of pulse in ms
        t2 : float
            end of pulse in ms. Only used for sinc and square pulses.
        OmegaM : float
            2*pi*(Rabi frequency) in kHz.
        type : str
            type of pulse to use. Current implement pulses are "square", "sinc" and "reburp"
        splitter : bool
            True for sinc in splitter configuration. False otherwise. Only used for "sinc".
        Return
        -----------------------
        Pulse value at time t
        """

        # transfrom t to numpy array
        t = np.array([t],dtype=float)
        
        # compute pulse
        if type == "square":
            return pulses.SquarePulse(t,t1,t2)[0]       
        elif type == "sinc":
            return pulses.SincPulse(t,t1,t2,OmegaM,splitter)[0]
        elif type == "reburp":
            return pulses.ReburpPulse(t,t1,OmegaM)[0]
        else:
            print("Waring: invalid pulse shape! Setting pulse equal to zero.")
            return 0.0
   
    def evolutionMatrix(self,tau,xi):
        """ Evolution matrix of the differential equation in dimensionless units 
        Parameters
        ---------------------
        tau : float
            dimensionless time at which matrix is computed
        xi : numpy array
            dimensionless momentum values at which matrix is to computed
        Return
        ------------------------
        Complex evolution matrix at time t.
        The matrix has shape (2*Nmax+1,2*Nmax+2,len(xi)) where Nmax = self.par["diffraction orders"]
        The orders as sorted as -Nmax,-Nmax+1,...,0,...,Nmax+1,Nmax. In other words index zero os -Nmax
        """

        """ compute useful quantities """

        Nmax = self.par["diffraction orders"] # number of orders to diffract
        nuR = self.par["Rabi frequency"] # Rabi frequency
        detuning = self.par["Detuning"] # detuning 
        slope = self.par["Phase slope"] # slope of frequency sweep
        t1 = self.par["Pulse beggining"] # beggining of pulse
        t2 = self.par["Pulse beggining"] + self.par["Pulse duration"] # end of pulse
        kB = self.par["Bragg wavevector"]
        g = self.par["gravity"]
        
        orders = np.arange(-Nmax,Nmax+1) # array with all orders
        OmegaM = 2*np.pi*nuR
        # compute pulse
        pulse = self.gpacket(tau/nuR,t1,t2,OmegaM,self.par["Pulse type"],self.par["Splitter"])

        # compute pulse phase and phase due to gravity
        phase = self.braggPhase(tau/nuR,detuning,slope) # pulse phase
        arg = phase + kB*g/2*np.power(tau/nuR,2)

        """ compute matrix elements """

        xiOnes = np.ones(len(xi),dtype=complex)

        # compute diagonal term, n = m
        A1 = np.diag(np.ones(2*Nmax+1,dtype=complex))
        A1 = np.tensordot(A1,xiOnes,axes = 0)

        # compute zero-th order frequency difference
        delta0 = self.par["Bragg recoil frequency"]/nuR*(2*xi-1)
        delta0 = np.tensordot(np.ones(len(orders),dtype=complex),delta0,axes=0)

        # compute term n-1 = m
        A20 = np.diag(np.ones(2*Nmax,dtype=complex),-1)*np.exp(-1j*arg)
        deltan = 2*(orders+1)*self.par["Bragg recoil frequency"]/nuR
        deltan = np.tensordot(deltan,xiOnes,axes=0)
        deltan = delta0 + deltan
        A21 = np.exp(1j*deltan*tau)
        A2 = np.einsum('nm,mk->nmk',A20,A21)

        # compute term n+1 = m
        A30 = np.diag(np.ones(2*Nmax,dtype=complex),1)*np.exp(1j*arg)
        deltan = 2*orders*self.par["Bragg recoil frequency"]/nuR
        deltan = np.tensordot(deltan,xiOnes,axes=0)
        deltan = delta0 + deltan
        A31 = np.exp(-1j*deltan*tau)
        A3 = np.einsum('nm,mk->nmk',A30,A31)

        return 2*np.pi*pulse*A1 + np.pi*pulse*(A2 + A3)
    
    def auxFunc(self,tau,xi,psi):
        """ Auxiliar function for self.braggSolver(). Basically here we perform the matrix multiplication 
        Parameters
        ---------------------
        tau : float
            dimensionless time at equation is to be solved
        xi : numpy array
            dimensionless momentum values at which matrix is to computed
        psi : wavefunction at time tau. psi should have dimension (2*Nmax+1,len(xi)) 
        where Nmax = self.par["diffraction orders"]
        Return
        ------------------------
        derivative of psi at time tau. It will have dimension (2*Nmax+1,len(xi))
        """

        # compute evolution matrix
        A = self.evolutionMatrix(tau,xi)

        # compute matrix multiplication
        return -1j*np.einsum('nmk,mk->nk',A,psi)
    
    def braggSolver(self,tau,xi,psi0):
        """ Solve the differential equation for Bragg diffraction 
        Parameters
        ---------------------
        tau : numpy array
            dimensionless time at equation is to be solved
        xi : numpy array
            dimensionless momentum values at which matrix is to computed
        psi0 : wavefunction at inital time tau. psi should have dimension (2*Nmax+1,len(xi)) 
        where Nmax = self.par["diffraction orders"]
        the order of psi0 should be -Nmax,-Nmax+1,...,0,...,Nmax+1,Nmax
        Return
        ------------------------
        psi at time tau. psi will have dimension (len(tau),2*Nmax+1,len(xi))
        the order of psi is -Nmax,-Nmax+1,...,0,...,Nmax+1,Nmax. In other words index zero is -Nmax
        """

        # number of orders to diffract
        Nmax = self.par["diffraction orders"] 
        # time step to solve differential equation
        step = self.par["time step solver"]*self.par["Rabi frequency"] 

        # create wavefunction
        psi = np.zeros((len(tau),2*Nmax+1,len(xi)),dtype=complex)
        
        # add initial time
        psi[0] = psi0

        # for each time
        for j in tqdm(range(1,len(tau))):

            # elements evaluation
            psi1 = step*self.auxFunc(tau[j-1],        xi, psi[j-1])
            psi2 = step*self.auxFunc(tau[j-1]+step/2, xi, psi[j-1]+psi1/2)
            psi3 = step*self.auxFunc(tau[j-1]+step/2, xi, psi[j-1]+psi2/2)
            psi4 = step*self.auxFunc(tau[j-1]+step,   xi, psi[j-1]+psi3)

            psi[j] = psi[j-1] + (psi1 + 2*psi2 + 2*psi3 + psi4)/6

        return psi    
    
    def compute_psiVelocity(self,t0,tfinal,v,psi0):
        """ Computes all the orders of the wavefunction in velocity space.
        Parameters
        ---------------------
        t0 : float
            initial time in ms
        tfinal : float
            final time in ms 
        v : numpy array
            velocity in mm/ms at which function is computed.
        psi0 : normalized wavefunction at inital time t0. psi should have dimension (2*Nmax+1,len(v)) 
        where Nmax = self.par["diffraction orders"]
        Return
        ------------------------
        t : array
            time values
        psi : array
            wavefunction at time t. psi will have dimension (len(t),2*Nmax+1,len(v)) where Nmax = self.par["diffraction orders"] 
            the order of psi is -Nmax,-Nmax+1,...,0,...,Nmax+1,Nmax. In other words index zero os -Nmax
        """

        # update parameters
        self.update_parameters()

        """ normalize quantities """

        # get times of pulse
        t1 = self.par["Pulse beggining"] # beggining of the pulse 
        t2 = t1 + self.par["Pulse duration"] # end of pulse

        # create time array
        t = np.arange(t0,t1,self.par["time step propagator"])
        t = np.concatenate((t,np.arange(t1,t2+self.par["time step solver"],self.par["time step solver"])))
        t = np.concatenate((t,np.arange(t2+self.par["time step solver"],tfinal,self.par["time step propagator"])))

        # transform time to dimensionless units
        tau = t*self.par["Rabi frequency"]
        # transform velocity to dimensionless units
        xi = self.par["mHe/hbar"]*v/self.par["Bragg wavevector"]

        """ initialize wavefunction """
        # number of orders to diffract
        Nmax = self.par["diffraction orders"] 
        # initialize psi
        psi = np.zeros((len(tau),2*Nmax+1,len(v)),dtype=complex)
        psi[0] = psi0

        """ if pulse is not on """
        mask = (t > t0)*(t < t1)
        psi[mask] = psi0

        """ if pulse is on we use self.braggSolver() to compute wavefunction """
        mask = (t >= t1)*(t<= t2)
        # compute wavefunction when pulse is on
        psi[mask] = self.braggSolver(tau[mask],xi,psi0)
        psi0 = psi[mask][-1]

        """ if pulse is not on """
        mask = (t > t2)*(t<=tfinal)
        psi[mask] = psi0

        return t , psi

    def createInitalPsi(self,v,psi0,m):
        """
        Given the inital wavefunction psi0 as a function of v, it creates the 2D array with dimension (2*Nmax+1,len(v)) 
        where Nmax = self.par["diffraction orders"] for the self.compute_psiVelocity()
        Parameters:
        ---------------------------------------
        v : 1D-numpy array
            array with velocity values
        psi0 : 1D-numpy array
            array with inital wavefunction as function of v at order m.
        m : int
            inital populated order. m should be between [-Nmax,Nmax].
        """
        # number of orders to diffract
        Nmax = self.par["diffraction orders"] 
        # make sure m is an int
        m = int(m)
        # if user chose an invalid order
        if np.abs(m) > Nmax:
            print("Invalid order. m must be smaller or equal to number of diffracted orders.")
            return 0 
        # create 2D array with dimension (2*Nmax+1,len(v)) 
        psi = np.zeros((2*Nmax+1,len(v)),dtype=complex)
        # populate order m with inital wavefunction
        psi[Nmax + m] = psi0
        return psi

    def nthPsi(self,psi,m):
        """ Given a psi ND-array with dimension (len(t),2*Nmax+1,...), or a 2D-array with dimension (len(t),2*Nmax+1)
        gets the m_th diffraction order wavefunction. v could be velocity or position, for example.
        Parameters
        ---------------------
        psi : ND-array or 2D-array
            It should have dimensions:
                (len(t),2*Nmax+1,...) or (len(t),2*Nmax+1)
            where Nmax = self.par["diffraction orders"] 
            The order of psi is -Nmax,-Nmax+1,...,0,...,Nmax+1,Nmax. In other words index zero os -Nmax
        m : float
            order of diffraction to get
        Return
        ------------------------
        psi : 2D-array or 1D-array
            Diffraction order m with dimensions
            (len(t),len(v)) or (len(t)) 
            where Nmax = self.par["diffraction orders"] 
            The order of psi is -Nmax,-Nmax+1,...,0,...,Nmax+1,Nmax. In other words index zero os -Nmax
        """
        # number of orders to diffract
        Nmax = self.par["diffraction orders"] 
        # make sure m is an int
        m = int(m)
        # if user chose an invalid order
        if np.abs(m) > Nmax:
            print("Invalid order. m must be smaller or equal to number of diffracted orders.")
            return 0 
        
        return psi[:,Nmax+m]
    
    def compute_psiSpace(self,psi,time,velocity,position,CP):
        """
        Computes the wavefunction in position space, given a psi in velocity space with 
        dimensions (len(time),2*Nmax+1,len(velocity)) where Nmax = self.par["diffraction orders"], 
        such that the order of psi is -Nmax,-Nmax+1,...,0,...,Nmax+1,Nmax. 
        This method is also vectorized, but depending on the size psi and position arrays, the calculation may require a lot RAM and
        crash the computer. Thus, you can split the calculation into smaller tensor multiplications using the variable CP.
        Parameters
        ---------------------
        psi : 3D-array
            It should have dimensions:
                (len(t),2*Nmax+1,len(velocity)) 
            where Nmax = self.par["diffraction orders"] 
            The order of psi is -Nmax,-Nmax+1,...,0,...,Nmax+1,Nmax. In other words index zero os -Nmax
        time : 1D numpy array
            array with time values in ms
        velocity : 1D numpy array
            array with velocity values in mm/ms
        position : 1D numpy array
            array with position values in mm at which wavefunction is to be determined
        CP : dict
            dictionary with two entris: CP["time"] and CP["position"]. Each entry tells how many times we split the time and position arrays
            to ligthen the calculation. 
        Return
        ------------------------
        psiSpace : 3D-array
            wavefunction in position space with dimensions:
                (len(t),2*Nmax+1,len(position))
            where Nmax = self.par["diffraction orders"] 
            The order of psi is -Nmax,-Nmax+1,...,0,...,Nmax+1,Nmax. In other words index zero os -Nmax
        """

        """ usefull parameters """
        # Bragg wavector
        kb = self.par["Bragg wavevector"]
        # rabi frequency
        nuR = self.par["Rabi frequency"]
        # number of orders to diffract
        Nmax = self.par["diffraction orders"] 
        # Bragg frequency
        omegab = self.par["Bragg recoil frequency"]
        # gravity
        g = self.par["gravity"]

        """ compute quantities in dimensionless units  """

        # compute dimensionless position
        x = position*kb/(2*np.pi)
        # compute dimensionless time
        tau = time*nuR
        # transform velocity to dimensionless units
        xi = self.par["mHe/hbar"]*velocity/kb
        # orders of diffraction
        orders = np.arange(-Nmax,Nmax+1)

        # normalize wavefunction to dimensionless units
        norm = np.sqrt(kb/self.par["mHe/hbar"])
        psi = psi*norm

        """ Compute function in position space """

        # split tau and x arrays into CP smaller arrays
        tauLists = np.array_split(tau,CP["time"])
        xLists = np.array_split(x,CP["position"])

        # initialize wavefucntion
        psiSpace = np.zeros((len(tau),2*Nmax+1,len(x)),dtype = complex)

        # for each tau subarray
        min = 0
        for a in tqdm(range(0,len(tauLists))):
            # for each x subarray
            xMin = 0
            for p in range(0,len(xLists)):

                """ transform variable arrays to shape (len(time),2*Nmax+1,len(velocity),len(position))"""

                # initalize base array for computation
                base = np.ones((len(tauLists[a]),2*Nmax+1,len(xi),len(xLists[p])),dtype=complex)
                # transfrom tau
                tau_N = np.einsum("ijkl,i->ijkl",base,tauLists[a])
                # transform orders
                orders_N = np.einsum("ijkl,j->ijkl",base,orders)
                # transfrom xi
                xi_N = np.einsum("ijkl,k->ijkl",base,xi)
                # transfrom x
                x_N = np.einsum("ijkl,l->ijkl",base,xLists[p])
                # transfrom psi
                psi_N = np.einsum("ijkl,ijk->ijkl",base,psi[min:min+len(tauLists[a])])

                """ commpute terms of integral """

                # compute phase
                phase = -omegab/nuR*np.power(xi_N+orders_N,2)*tau_N
                phase = phase + (g*kb/(2*np.power(nuR,2))*np.power(tau_N,2)+2*np.pi*x_N)*(xi_N+orders_N)
                # compute argument of integral
                psi_N = psi_N*np.exp(1j*phase)

                """ perform integration over xi """

                xi = self.par["mHe/hbar"]*velocity/kb
                psiSpace[min:min+len(tauLists[a]),:,xMin:xMin+len(xLists[p])] = integrate.simpson(y = psi_N , x = xi , axis = 2)
    
                # update index
                xMin = xMin + len(xLists[p])

            # update index
            min = min + len(tauLists[a])

        return psiSpace*np.sqrt(kb/(2*np.pi))
    
    
    # delete after testing finished
    def test(self,tau,xi):

        Nmax = self.par["diffraction orders"] # number of orders to diffract
        nuR = self.par["Rabi frequency"] # Rabi frequency
        detuning = self.par["Detuning"] # detuning 
        slope = self.par["Phase slope"] # slope of frequency sweep
        t1 = self.par["Pulse beggining"] # beggining of pulse
        t2 = self.par["Pulse beggining"] + self.par["Pulse duration"] # end of pulse
        kB = self.par["Bragg wavevector"]
        g = self.par["gravity"]
        
        orders = np.arange(-Nmax,Nmax+1) # array with all orders
        OmegaM = 2*np.pi*nuR
        # compute pulse
        pulse = self.Gpacket(tau/nuR,t1,t2,OmegaM,self.par["Pulse type"],self.par["Splitter"])

        # compute pulse phase and phase due to gravity
        phase = self.BraggPhase(tau/nuR,detuning,slope) # pulse phase
        arg = phase + kB*g/2*np.power(tau/nuR,2)

        A = np.zeros((2*Nmax+1,2*Nmax+1,len(xi)),dtype = complex)

        """ compute matrix elements """
        for n in range(-Nmax,Nmax+1,1):
            for m in range(-Nmax,Nmax+1,1):
                for k in range(0,len(xi)):
                    delta0 = self.par["Bragg recoil frequency"]/nuR*(2*xi[k]-1)
                    # for diagonal terms
                    if n == m:
                        A[n+Nmax,m+Nmax,k] = 2*np.pi*pulse
                        pass
                    elif (n-1) == m:
                        deltan =  delta0 + 2*(m+1)*self.par["Bragg recoil frequency"]/nuR
                        A[n+Nmax,m+Nmax,k] = np.pi*pulse*np.exp(-1j*(arg-deltan*tau))
                    elif (n+1) == m:
                        deltan =  delta0 + 2*m*self.par["Bragg recoil frequency"]/nuR
                        A[n+Nmax,m+Nmax,k] = np.pi*pulse*np.exp(1j*(arg-deltan*tau))
                        pass
                    else:
                        pass

        return A








        




        
        
