# ----------------------------------
# Created on the 09/2025 by Rui
#
# Copyright (c) 2025 - Helium1@LCF
# ----------------------------------
#

import sys
import numpy as np
import scipy.integrate as integrate
from scipy.interpolate import RegularGridInterpolator
from tqdm import tqdm

""" Units of quantities with dimensions are mm/ms/kg.

We use Runge-kutta 4th order to solve the differential equation. """

class Bragg:
    
    """   Object initialization  """
    def __init__(self,vB):
        """ vB is the Bragg velocity in mm/ms. The class has a dictionary self.par with all relevant parameters 
        to the experiment and calculation. In case you modify "element A" of self.par, you 
        can update the elements of self.par that depend on "element A" using the method self.update_parameters() """
        
        # dictionary to hold all relevant parameters in mm/ms/kg
        self.par = dict()

        """ Define constants in mm/ms/kg units """
        self.par["hbar"] = 1.054*10**(-31) # mm^2 kg/ms
        self.par["mHe"] = 6.65*10**(-27) # Helium's mass in kg
        self.par["gravity"] = 9.8067e-3 # mm/ms^2

        """ Define experimental lattice parameters in mm/ms/kg units """
        self.par["Bragg velocity"] = vB # Bragg velocity
        self.par["Rabi frequency"] = 1 # Rabi frequency
        self.par["Detuning"] = 10.2 # define detuning in kHz
        self.par["Pulse beggining"] = 0.0 # beggining of the pulse 
        self.par["Pulse duration"] = 1 # duration of the pulse
        self.update_parameters() # update dictionary
        
        """ Define calculation parameters """
        self.par["velocity step"] = 1e-4 # define velocity step in mm/ms
        self.par["diffraction orders"] = 5 # define number of diffraction orders
        # define propagation step when pulse is not on
        self.par["time step propagator"] = 0.1
        # define time step to solve differential equation in ms when pulse is on
        if self.par["Rabi frequency"] > self.par["Bragg recoil frequency"]/(2*np.pi):
            self.par["time step solver"] = 1/(20*self.par["Rabi frequency"])
        else:
            self.par["time step solver"] = 1/(20*self.par["Bragg recoil frequency"]/(2*np.pi))

        """ atribute to turn on or off tqdm bar """
        self.tqdm = False
        
        self.update_parameters() # update dictionary
        
    def update_parameters(self):
        """ Updates self.par which is a dictionary with all relevant experimental and calculation parameters """
        
        # update self.par
        self.par["mHe/hbar"] = self.par["mHe"]/self.par["hbar"]
        self.par["Bragg wavevector"] = self.par["mHe/hbar"]*self.par["Bragg velocity"]
        self.par["Bragg recoil frequency"] = np.power(self.par["Bragg wavevector"],2)/(2*self.par["mHe/hbar"])
        self.par["Frequency slope to compensate gravity"] = -self.par["Bragg wavevector"]*self.par["gravity"]/(2*np.pi)

    def frequencyOrderDifference(self,n,v):
        """ Frequency order difference (2 pi factor included) \delta_n = (\epsilon_{n+1} - \epsilon_{n})/hbar
        Parameters
        ------------------------
        n : float
            diffraction order to consider. E.g if n = 1, than FrequencyOrderDifference() gives the frequency difference between
            order n = 1 and order n = 0.
        v : float
            velocity of the atom in mm/ms units
        Return 
        -------------------------------------------
            The 2*np.pi*frequency difference between diffraction orders
        """
        # update parameter
        self.update_parameters()

        delta0 = v*self.par["Bragg wavevector"] - self.par["Bragg recoil frequency"]
        return delta0 + 2*n*self.par["Bragg recoil frequency"]
    
    def initializePulse(self,t,pulse,phase):
        """ Function that interpolates the user defined pulse and phase for calculation. After interpolation, we have two functions,
        self.gpacket(t) and self.braggPhase(t) used for calculation
        Parameters
        -----------------------------------
        t : numpy array
            time in ms
        pulse : numpy array
            pulse values at time t
        phase : numpy array
            phase of the pulse at time t
        """

        self.gpacket = RegularGridInterpolator([t],pulse,method = "cubic")
        self.braggPhase = RegularGridInterpolator([t],phase,method = "cubic")
   
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
        kB = self.par["Bragg wavevector"]
        g = self.par["gravity"]
        
        orders = np.arange(-Nmax,Nmax+1) # array with all orders
        # compute pulse
        pulse = self.gpacket([tau/nuR])[0]

        # compute pulse phase and phase due to gravity
        phase = self.braggPhase([tau/nuR])[0] # pulse phase
        arg = phase + kB*g/2*np.power(tau/nuR,2) # add phase due to gravity

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
    
    def braggSolver(self,tau_i,tau_f,xi,psi0,saveEvery):
        """ Solve the differential equation for Bragg diffraction 
        Parameters
        ---------------------
        tau_0 : float
            dimensionless inital time at equation is to be solved
        tau_f : float
            dimensionless final time
        xi : numpy array
            dimensionless momentum values at which matrix is to computed
        psi0 : 2D-numpy array
            wavefunction at inital time tau. psi should have dimension (2*Nmax+1,len(xi)) 
            where Nmax = self.par["diffraction orders"]
            the order of psi0 should be -Nmax,-Nmax+1,...,0,...,Nmax+1,Nmax
        saveEvery : float
            save calculation result of differential equation every saveEvery dimensionless time
        Return
        ------------------------
        time : 1D-numpy array
            array with time values
        psi : 3D-numpy complex array
            psi at time time. psi will have dimension (len(time),2*Nmax+1,len(xi))
        the order of psi is -Nmax,-Nmax+1,...,0,...,Nmax+1,Nmax. In other words index zero is -Nmax
        """

        # time step to solve differential equation
        step = self.par["time step solver"]*self.par["Rabi frequency"] 
        # create time array for computation
        tau = np.arange(tau_i,tau_f+step,step)
        # initialize time array to save data
        time = []
        # initialize wavefunction
        psi = []
        
        # add initial time
        time.append(tau[0])
        psi.append(psi0)

        # initialize timer
        timer = 0.0
        # for each time
        for j in tqdm(range(1,len(tau)),disable = self.tqdm):

            # increment timer
            timer = timer + step

            # elements evaluation
            psi1 = step*self.auxFunc(tau[j-1],        xi, psi0)
            psi2 = step*self.auxFunc(tau[j-1]+step/2, xi, psi0+psi1/2)
            psi3 = step*self.auxFunc(tau[j-1]+step/2, xi, psi0+psi2/2)
            psi4 = step*self.auxFunc(tau[j-1]+step,   xi, psi0+psi3)

            # compute wavefunction at time tau[j]
            psi0 = psi0 + (psi1 + 2*psi2 + 2*psi3 + psi4)/6

            # save every saveEvery
            if timer >= saveEvery:
                # save values
                time.append(tau[j])
                psi.append(psi0)
                # reset timer
                timer = 0.0

        # append last result of calculation
        time.append(tau[-1])
        psi.append(psi0)

        return np.array(time,dtype = float)/self.par["Rabi frequency"] , np.array(psi,dtype=complex)  
    
    def compute_psiVelocity(self,t0,tfinal,v,psi0,saveEvery):
        """ Computes all the orders of the wavefunction in velocity space.
        Parameters
        ---------------------
        t0 : float
            initial time in ms
        tfinal : float
            final time in ms 
        v : numpy array
            velocity in mm/ms at which function is computed.
        psi0 : 2D-numpy array
            normalized wavefunction at inital time t0. psi should have dimension (2*Nmax+1,len(v)) 
            where Nmax = self.par["diffraction orders"]
        saveEvery : float
            save calculation result of differential equation every saveEvery ms
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

        # number of orders to diffract
        Nmax = self.par["diffraction orders"] 

        """ if pulse is not on """
        # create time array
        time_1 = np.arange(t0,t1,self.par["time step propagator"])
        # create wavefunction
        psi_1 = np.zeros((len(time_1),2*Nmax+1,len(v)),dtype=complex)
        # add initial condition
        psi_1[:] = psi0

        """ if pulse is on we use self.braggSolver() to compute wavefunction """
        # transform velocity to dimensionless units
        xi = self.par["mHe/hbar"]*v/self.par["Bragg wavevector"]
        # compute intial and final times in dimensionless units
        tau_i = t1*self.par["Rabi frequency"]
        tau_f = t2*self.par["Rabi frequency"]
        # compute wavefunction when pulse is on
        time_2 , psi_2 = self.braggSolver(tau_i,tau_f,xi,psi0,saveEvery*self.par["Rabi frequency"])

        """ if pulse is not on """
        # create time array
        time_3 = np.arange(time_2[-1],tfinal+self.par["time step propagator"],self.par["time step propagator"])
        # create wavefunction
        psi_3 = np.zeros((len(time_3),2*Nmax+1,len(v)),dtype=complex)
        # add final condition
        psi_3[:] = psi_2[-1]

        """ concatenate results """
        psi = np.concatenate((psi_1,psi_2,psi_3),axis = 0)
        del psi_1
        del psi_2
        del psi_3
        t = np.concatenate((time_1,time_2,time_3))
        del time_1
        del time_2
        del time_3

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
        Return
        ------------------------------------
            Inital wavefunction
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
        m : int
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
        """ Computes the wavefunction in position space, given a psi in velocity space with 
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
        for a in tqdm(range(0,len(tauLists)),disable = self.tqdm):
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
                phase = -omegab/nuR*tau_N*np.power(xi_N+orders_N,2)
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
    
    def saveToFile(self,filename,psi,time,coordinates):
        """ Saves to file.npy a wavefunction.
        Parameters
        ---------------------------------------
        filename : str
            name of file
        psi : 3D-numpy array
            wavefunction to save.
        time : 1D-numpy array
            array with time values
        coordinates : 1D-numpy array
            coordinates values at which psi is known
        """

        with open(filename, 'wb') as f:
            np.save(f, psi)
            np.save(f, time)
            np.save(f, coordinates)

        print("Data saved to file!")
    
    def LoafFile(self,filename):
        """ Loads wavefunction from file.npy.
        Parameters
        ---------------------------------------
        filename : str
            name of file
        Return 
        ---------------------------------------
        psi : 3D-numpy array
            wavefunction to save.
        time : 1D-numpy array
            array with time values
        coordinates : 1D-numpy array
            coordinates values at which psi is known
        """

        with open(filename, 'rb') as f:
            psi = np.load(f)
            time = np.load(f)
            coordinates = np.load(f)

        return psi , time , coordinates

    def MemorySize(self,shape):
        """ Computes the memory size of a complex array 
        Parameters 
        ------------------------------------------------
        shape : ND-array or N-tuple
            shape of array to compute size
        Return 
        ------------------------------------------------
            size of array in Giga-bytes
        """
        # size of a complex number in bytes
        size = 16.0
        for i in shape:
            size = size*i
        return size/1073741824

        




        
        
