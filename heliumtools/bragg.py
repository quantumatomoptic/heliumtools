import numpy as np
import matplotlib.pyplot as plt
from odeintw import odeintw
import scipy.constants as cst
from tqdm import tqdm


hbar = cst.hbar
m = 4 * cst.physical_constants["atomic mass constant"][0]
pi = np.pi


class Pulse:
    def set_params(self,params):
        self.wavevector = params["wavevector"]
        self.amplitude = params["amplitude"]
        self.duration = params["duration"]
        self.type = params["resonance frequencies"]#1 freq or 2 freq
        self.waveform = params["shape"]
        self.pulse_type = params["pulse type"]
        self.window_coefficient = params["window coefficient"]
        self.window_delay = params["window delay"]
        self.phase_difference = params["phase difference"]
        self.phase = params["global phase"]
        self.detuning = params["detuning"]
        
    def window(self,t):
        T = self.duration
        T_delay = self.window_delay
        alpha = self.window_coefficient
        if t < T_delay or t > T-T_delay:       
            return 0
        else:
            return((np.sin(pi*(t-T_delay)/(T-2*T_delay)))**(alpha))          
        

    def modulation(self,t):
        if self.type == 1:
            return 1.0
        elif self.type == 2:
            OmegaD = self.detuning
            T = self.duration
            theta = self.phase_difference
            return np.exp(+1j*OmegaD*(t-T/2)/2+1j*theta/2)+np.exp(-1j*OmegaD*(t-T/2)/2-1j*theta/2)
        
    def shape(self,t):
        if self.waveform == "constant":
            return 1.0
        elif self.waveform == "sinc":
            Omega_M = self.amplitude 
            if self.pulse_type == "mirror":
                Omega_S = Omega_M
            elif self.pulse_type == "splitter":
                Omega_S = 2 * Omega_M
            T = self.duration
            return np.sin(Omega_S*(t-T/2))/(Omega_S*(t-T/2))
    
    def pulse(self,t):
        if t<0 or t>self.duration:
            return 0
        else:
            Omega_M = self.amplitude    
            return Omega_M * self.shape(t) * self.modulation(t) * self.window(t) 
        
    def plot(self):
        number_of_points = 500
        T = self.duration
        time = np.linspace( 0 - 0.1 * T , T * 1.1 , number_of_points)
        pulse_abs = []
        for k in range(number_of_points):
            pulse_abs.append(np.abs(self.pulse(time[k])))
            
        plt.figure(figsize=(10,7))
        plt.grid(True)
        plt.xlabel("Time (µs)")
        plt.ylabel(r"$|\Omega_R|$ (kHz)")
        plt.plot(time * 1e6, np.array(pulse_abs)/(2*pi)*1e-3)
        plt.show()

        

       
class Bragg():

    
# SET PARAMETERS
    def set_params(self,params,pulse):
        self.levels = params["number of levels"]
        self.detuning = params["detuning"]
        self.min_detuning = params["min detuning"]
        self.max_detuning = params["max detuning"]
        self.number_of_detuning_points = params["number of points"]
        self.boxes_number = params["number of boxes"]
        self.box_width = params["box width"]
        self.show_boxes = params["show boxes"]
        self.duration = pulse.duration
        self.pulse_type = pulse.pulse_type
        self.number_of_points = 1000
        T = self.duration
        self.time = np.linspace( 0, T , self.number_of_points)
        self.detuning_list = np.linspace(self.min_detuning, self.max_detuning, self.number_of_detuning_points)


    
# SOLVER
    def time_solver(self,Bragg_pulse):
        
        def pulse(t):
            return Bragg_pulse.pulse(t)
        
        def zfunc2(z, t):
            z1, z2 = z
            return [z2*1j*pulse(t)*np.exp(1j*delta*t)/2, z1*1j*np.conj(pulse(t))*np.exp(-1j*delta*t)/2]

        def zjac2(z, t):
            z1, z2 = z
            jac = np.array([[0, 1j*pulse(t)*np.exp(1j*delta*t)/2], [1j*np.conj(pulse(t))*np.exp(-1j*delta*t)/2, 0]])
            return jac
        
        #sets the detuning
        delta = self.detuning
        #initial state
        z0 = np.array([complex(1.0), complex(0.0)])
        #solve the equation
        z, infodict = odeintw(zfunc2, z0, self.time, Dfun=zjac2,full_output=True)
        #return the coefficients
        C0=z[:,0]
        C2=z[:,1]
        return(C0,C2)
  


    def time_solver_5lvl(self,Bragg_pulse):
        k = Bragg_pulse.wavevector #a Brugg pulse transfers 2 * hbar * k
        Omega_recoil = hbar * ( k**2 ) / ( 2 * m ) #1 photon recoil energy
        OmegaB = 4 * Omega_recoil

        def pulse(t):
            return Bragg_pulse.pulse(t)
        
        def zfunc2(z, t):
            z1, z2, z3, z4, z5 = z
            return [z2*1j*pulse(t)*np.exp(1j*(delta+4*OmegaB)*t)/2, z1*1j*np.conj(pulse(t))*np.exp(-1j*(delta+4*OmegaB)*t)/2+z3*1j*pulse(t)*np.exp(1j*(delta+2*OmegaB)*t)/2,z2*1j*np.conj(pulse(t))*np.exp(-1j*(delta+2*OmegaB)*t)/2+z4*1j*pulse(t)*np.exp(1j*(delta)*t)/2,z3*1j*np.conj(pulse(t))*np.exp(-1j*(delta)*t)/2+z5*1j*pulse(t)*np.exp(1j*(delta-2*OmegaB)*t)/2,z4*1j*np.conj(pulse(t))*np.exp(-1j*(delta-2*OmegaB)*t)/2]

        def zjac2(z, t):
            z1, z2, z3, z4, z5 = z
            jac = np.array([[0, 1j*pulse(t)*np.exp(1j*(delta+4*OmegaB)*t)/2,0,0,0], [1j*np.conj(pulse(t))*np.exp(-1j*(delta+4*OmegaB)*t)/2, 0,1j*pulse(t)*np.exp(1j*(delta+2*OmegaB)*t)/2,0,0],[0,1j*np.conj(pulse(t))*np.exp(-1j*(delta+2*OmegaB)*t)/2, 0,1j*pulse(t)*np.exp(1j*(delta)*t)/2,0],[0,0, 1j*np.conj(pulse(t))*np.exp(-1j*(delta)*t)/2,0,1j*pulse(t)*np.exp(1j*(delta-2*OmegaB)*t)/2],[0,0,0, 1j*np.conj(pulse(t))*np.exp(-1j*(delta-2*OmegaB)*t)/2,0]])
            return jac
        
        #sets the detuning
        delta = self.detuning
        #initial state
        z0 = np.array([complex(0.0),complex(0.0),complex(1.0), complex(0.0),complex(0.0)])
        #solve the equation
        z, infodict = odeintw(zfunc2, z0, self.time, Dfun=zjac2,full_output=True)
        #return the coefficients
        C0=z[:,2]
        C2=z[:,3]
        return(C0,C2)
            
        
    def compute(self,pulse):
        T = self.duration
        P_C0 = []
        P_C2 = []
        a_C0 = []
        a_C2 = []
        for k in tqdm(range(self.number_of_detuning_points), desc="Progress"):
        #for k in range(self.number_of_detuning_points):
            self.detuning = self.detuning_list[k]
            if self.levels == 2:
                C0, C2 = self.time_solver(pulse)
            elif self.levels == 5 :
                C0, C2 = self.time_solver_5lvl(pulse)
            P_C0.append((np.abs(C0)**2)[self.number_of_points-1])
            P_C2.append((np.abs(C2)**2)[self.number_of_points-1])
            a_C0.append(np.angle(C0[self.number_of_points-1]))
            a_C2.append(np.angle(C2[self.number_of_points-1]) + self.detuning_list[k] * self.duration / 2 % (2*pi))
        return P_C0, P_C2, a_C0, a_C2

    
# PLOT FUNCTIONS
    # Time plot
    def time_plot_population(self,pulse):
        C0, C2 = self.time_solver(pulse)
        plt.figure(figsize=(10,7))
        plt.grid(True)
        plt.xlabel("Time (ms)")
        plt.ylabel("Population") 
        plt.plot(self.time * 1e3, np.abs(C0)**2, label = r'$|C0|^2 $')
        plt.plot(self.time * 1e3, np.abs(C2)**2, label = r'$|C2|^2 $')
        plt.legend()
        plt.show()
        
    def time_plot_phase(self,pulse):
        C0, C2 = self.time_solver(pulse)
        plt.figure(figsize=(10,7))
        plt.grid(True)
        plt.xlabel("Time (ms)")
        plt.ylabel("Phase (rad)") 
        plt.plot(self.time * 1e3, np.angle(C0), label = 'arg(C0)')
        plt.plot(self.time * 1e3, np.angle(C2), label = 'arg(C2)')
        plt.legend()
        plt.show()  
        
    # Detuning plot
    def plot_boxes(self,plt):
        if self.show_boxes is True:
            for k in range(self.boxes_number+1):
                plt.axvline(x=-((2*k+1)/2)*self.box_width/(2*pi) *1e-3,color='black',linestyle='--',linewidth=1)
                plt.axvline(x=+((2*k+1)/2)*self.box_width/(2*pi) *1e-3,color='black',linestyle='--',linewidth=1)
            plt.axvspan(-self.box_width/(2*pi) *(1e-3)/2-self.boxes_number*self.box_width/(2*pi) *(1e-3), -self.box_width/(2*pi) *(1e-3)/2, alpha=0.2, color='red')
            plt.axvspan(self.box_width/(2*pi) *(1e-3)/2, self.box_width/(2*pi) *(1e-3)/2+self.boxes_number*self.box_width/(2*pi) *(1e-3),alpha=0.2, color='green')
    
    
    def plot_reflectivity(self,results):
        #P_C0, P_C2, a_C0, a_C2 = compute(pulse)
        plt.figure(figsize=(10,7))
        plt.grid(True)
        self.plot_boxes(plt)
        plt.xlabel("Detuning (kHz)")
        plt.ylabel("Population") 
        for k in range(len(results)):
            P_C0, P_C2, a_C0, a_C2, label_curve = results[k]
            plt.plot(self.detuning_list / (2*pi) *1e-3, P_C0, label = r'$|C_0|^2$'+label_curve)
            plt.plot(self.detuning_list / (2*pi) *1e-3, P_C2, label = r'$|C_2|^2$'+label_curve)
        plt.legend()
        plt.show()        

        
    
    def plot_phase(self,results):
        #P_C0, P_C2, a_C0, a_C2 = results
        plt.figure(figsize=(10,7))
        plt.grid(True)
        self.plot_boxes(plt)
        plt.xlabel("Detuning (kHz)")
        plt.ylabel("Angle (°)") 
        for k in range(len(results)):
            P_C0, P_C2, a_C0, a_C2, label_curve = results[k]
            reflected_angle = ((np.unwrap(4 * np.array(a_C2), discont = pi) / 4)) * 180 / pi
            transmitted_angle = ((np.unwrap(4 * np.array(a_C0), discont = pi) / 4)) * 180 / pi
            if self.pulse_type == "splitter":
                imprinted_phase = transmitted_angle - reflected_angle
            elif self.pulse_type == "mirror":
                imprinted_phase = - reflected_angle
            plt.plot(self.detuning_list / (2*pi) *1e-3, imprinted_phase, label = label_curve)
            plt.legend()
        plt.show()               
        
# TOOLS       
    def display(self):
        print(np.min(self.detuning_list / (2*pi) *1e-3))
        print(np.max(self.detuning_list / (2*pi) *1e-3))
        
        
        
        
