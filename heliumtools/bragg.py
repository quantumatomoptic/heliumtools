import numpy as np
import matplotlib.pyplot as plt
from odeintw import odeintw
import scipy.constants as cst
from tqdm import tqdm


#constants
hbar = cst.hbar
m = 4 * cst.physical_constants["atomic mass constant"][0]
pi = np.pi


#params
#plt.rcParams.update({'font.size': 14})

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
        self.drift = params["phase drift"]
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

    def phase_drift(self,t):
        drift = self.drift #value (in degrees) of the phase drift during the duration of the pulse
        omega_drift = drift / self.duration
        return np.exp( 1j * omega_drift * t )

    def phase_noise(self,t):
        drift = 2.5 * pi / 180 #value (in degrees) of the phase drift during the duration of the pulse
        omega_drift = 2 * pi * 3 / self.duration
        return np.exp( 1j * 2 * drift * np.sin( omega_drift * t ) )

    def global_phase(self,t):
        phase = 0 * pi / 180
        return np.exp( 1j * phase )

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
            return Omega_M * self.shape(t) * self.modulation(t) * self.window(t) * self.phase_drift(t) * self.global_phase(t)

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



class Bell():


    def set_params(self, params_Bell, pulse_1, pulse_2):
        self.mirror_center = params_Bell["mirror center"]
        self.splitter_center = params_Bell["splitter center"]
        self.mirror = pulse_1
        self.splitter = pulse_2
        #durations
        self.mirror_duration = pulse_1.duration
        self.splitter_duration = pulse_2.duration
        self.mirror_delay = self.mirror_center  - pulse_1.duration/2
        self.splitter_delay = self.splitter_center - pulse_2.duration/2 - (self.mirror_delay + pulse_1.duration)
        self.total_duration = self.mirror_delay +  self.mirror_duration + self.splitter_delay + self.splitter_duration
        #absolute times
        self.t_start_mirror = self.mirror_delay
        self.t_stop_mirror = self.mirror_delay + self.mirror_duration
        self.t_start_splitter = self.mirror_delay + self.mirror_duration + self.splitter_delay
        self.time = np.linspace(0,self.total_duration,1000)
        #detuning
        self.min_detuning = params_Bell["min detuning"]
        self.max_detuning = params_Bell["max detuning"]
        self.number_of_detuning_points = params_Bell["number of points"]
        self.detuning_list = np.linspace(self.min_detuning, self.max_detuning, self.number_of_detuning_points)
        self.boxes_number = params_Bell["number of boxes"]
        self.box_width = params_Bell["box width"]
        self.show_boxes = params_Bell["show boxes"]
        #parameters
        self.omega_0 = 0
        self.omega_2 = 0
        self.Delta_omega = 2 * pi * 25e3 #ie 50 mm/s
        self.remove_transmitted_atoms = params_Bell["remove transmitted atoms"]


    def Bell_pulse(self,t):
        if t < self.t_start_mirror :
            return 0
        elif t >= self.t_start_mirror and t < self.t_stop_mirror :
            return self.mirror.pulse( t - self.t_start_mirror )
        elif t >= self.t_stop_mirror and t < self.t_start_splitter :
            return 0
        elif t >= self.t_start_splitter :
            return self.splitter.pulse( t - self.t_start_splitter)


    def time_solver(self, initial_state, time):

        def pulse(t):
            return self.Bell_pulse(t)

        def zfunc2(z, t):
            z1, z2 = z
            return [-1j*omega0*z1+z2*1j*pulse(t)*np.exp(1j*(delta+Delta_omega)*t)/2, z1*1j*np.conj(pulse(t))*np.exp(-1j*(delta+Delta_omega)*t)/2-1j*omega2*z2]

        def zjac2(z, t):
            z1, z2 = z
            jac = np.array([[-1j*omega0, 1j*pulse(t)*np.exp(1j*(delta+Delta_omega)*t)/2], [1j*np.conj(pulse(t))*np.exp(-1j*(delta+Delta_omega)*t)/2, -1j*omega2]])
            return jac

        #sets the detuning
        delta = self.detuning
        omega0 = self.omega_0
        omega2 = self.omega_2
        Delta_omega = self.Delta_omega
        #initial state
        z0 = initial_state
        #solve the equation
        z, infodict = odeintw(zfunc2, z0, time, Dfun=zjac2,full_output=True)
        #return the coefficients
        C0=z[:,0]
        C2=z[:,1]
        return(C0,C2)


    def compute(self):

        #initial state |psi_in>= (|Alice1> x |Bob1> + e^{i*phi0} |Alice2> x |Bob2>)/sqrt(2)
        phi0=0
        Alice1_init=[complex(1.0),complex(0.0)]
        Bob1_init=[complex(0.0),complex(1.0)]
        Alice2_init=[complex(0.0),complex(1.0)]
        Bob2_init=[complex(1.0),complex(0.0)]


        #each list will contain the value of the final state for a given detuning
        Alice1_res=[]
        Alice2_res=[]
        Bob1_res=[]
        Bob2_res=[]

        if not self.remove_transmitted_atoms:
            for k in tqdm(range(self.number_of_detuning_points), desc="Progress"):
                # Alice doublet
                self.detuning = - self.detuning_list[k]
                self.omega_0 = - self.detuning_list[k]
                self.omega_2 = - self.detuning_list[k] + self.Delta_omega
                C0_Alice1, C2_Alice1 = self.time_solver(Alice1_init, self.time)
                Alice1_res.append([C0_Alice1[len(C0_Alice1)-1],C2_Alice1[len(C0_Alice1)-1]])
                C0_Alice2, C2_Alice2 = self.time_solver(Alice2_init, self.time)
                Alice2_res.append([C0_Alice2[len(C0_Alice2)-1],C2_Alice2[len(C0_Alice2)-1]])
                # Bob doublet
                self.detuning = + self.detuning_list[k]
                self.omega_0 = self.detuning_list[k]
                self.omega_2 = self.detuning_list[k] + self.Delta_omega
                C0_Bob1, C2_Bob1 = self.time_solver(Bob1_init, self.time)
                Bob1_res.append([C0_Bob1[len(C0_Bob1)-1],C2_Bob1[len(C0_Bob1)-1]])
                C0_Bob2, C2_Bob2 = self.time_solver(Bob2_init, self.time)
                Bob2_res.append([C0_Bob2[len(C0_Bob2)-1],C2_Bob2[len(C0_Bob2)-1]])

        else:
            time1 = np.linspace(0,self.t_start_mirror,300)
            time2 = np.linspace(self.t_start_mirror,self.t_stop_mirror,300)
            time3 = np.linspace(self.t_stop_mirror,self.total_duration,300)

            for k in tqdm(range(self.number_of_detuning_points), desc="Progress"):


                # Alice doublet
                self.detuning = - self.detuning_list[k]
                self.omega_0 = - self.detuning_list[k]
                self.omega_2 = - self.detuning_list[k] + self.Delta_omega
                # STEP 1 : FREE FALL 1
                C0_Alice1, C2_Alice1 = self.time_solver(Alice1_init, time1)
                lastval = len(C0_Alice1)-1
                Alice1 = [C0_Alice1[lastval],C2_Alice1[lastval]]
                C0_Alice2, C2_Alice2 = self.time_solver(Alice2_init, time1)
                Alice2 = [C0_Alice2[lastval],C2_Alice2[lastval]]
                # STEP 2 : MIRROR
                Alice1_a = [Alice1[0],0]
                _, C2_Alice1 = self.time_solver(Alice1_a, time2)
                Alice1_b = [0,Alice1[1]]
                C0_Alice1, _ = self.time_solver(Alice1_b, time2)
                lastval = len(C0_Alice1)-1
                Alice1 = [C0_Alice1[lastval],C2_Alice1[lastval]]
                Alice2_a = [Alice2[0],0]
                _, C2_Alice2 = self.time_solver(Alice2_a, time2)
                Alice2_b = [0,Alice2[1]]
                C0_Alice2, _ = self.time_solver(Alice2_b, time2)
                Alice2 = [C0_Alice2[lastval],C2_Alice2[lastval]]
                # STEP 3 : FREE FALL 2 + SPLITTER
                C0_Alice1, C2_Alice1 = self.time_solver(Alice1, time3)
                lastval = len(C0_Alice1)-1
                Alice1_res.append([C0_Alice1[lastval],C2_Alice1[lastval]])
                C0_Alice2, C2_Alice2 = self.time_solver(Alice2, time3)
                Alice2_res.append([C0_Alice2[lastval],C2_Alice2[lastval]])

                # Bob doublet
                # STEP 1 : FREE FALL 1
                self.detuning = + self.detuning_list[k]
                self.omega_0 = self.detuning_list[k]
                self.omega_2 = self.detuning_list[k] + self.Delta_omega
                C0_Bob1, C2_Bob1 = self.time_solver(Bob1_init, time1)
                Bob1 = [C0_Bob1[lastval],C2_Bob1[lastval]]
                C0_Bob2, C2_Bob2 = self.time_solver(Bob2_init, time1)
                Bob2 = [C0_Bob2[lastval],C2_Bob2[lastval]]
                # STEP 2 : MIRROR
                Bob1_a = [Bob1[0],0]
                _, C2_Bob1 = self.time_solver(Bob1_a, time2)
                Bob1_b = [0,Bob1[1]]
                C0_Bob1, _ = self.time_solver(Bob1_b, time2)
                Bob1 = [C0_Bob1[lastval],C2_Bob1[lastval]]
                Bob2_a = [Bob2[0],0]
                _, C2_Bob2 = self.time_solver(Bob2_a, time2)
                Bob2_b = [0,Bob2[1]]
                C0_Bob2, _ = self.time_solver(Bob2_b, time2)
                Bob2 = [C0_Bob2[lastval],C2_Bob2[lastval]]
                # STEP 3 : FREE FALL 2 + SPLITTER
                C0_Bob1, C2_Bob1 = self.time_solver(Bob1, time3)
                Bob1_res.append([C0_Bob1[lastval],C2_Bob1[lastval]])
                C0_Bob2, C2_Bob2 = self.time_solver(Bob2, time3)
                Bob2_res.append([C0_Bob2[lastval],C2_Bob2[lastval]])



        # final state can be written |psi_out>=A|-->+B|-+>+C|+->+D|++>
        A=[]
        B=[]
        C=[]
        D=[]

        #compute probabilities and correlator
        for k in range(len(self.detuning_list)):
            A.append((1/np.sqrt(2))*(Alice1_res[k][0]*Bob1_res[k][0]+np.exp(1j*phi0)*Alice2_res[k][0]*Bob2_res[k][0]))
            B.append((1/np.sqrt(2))*(Alice1_res[k][0]*Bob1_res[k][1]+np.exp(1j*phi0)*Alice2_res[k][0]*Bob2_res[k][1]))
            C.append((1/np.sqrt(2))*(Alice1_res[k][1]*Bob1_res[k][0]+np.exp(1j*phi0)*Alice2_res[k][1]*Bob2_res[k][0]))
            D.append((1/np.sqrt(2))*(Alice1_res[k][1]*Bob1_res[k][1]+np.exp(1j*phi0)*Alice2_res[k][1]*Bob2_res[k][1]))

        Bell_correlator=[]
        Ppp=[]
        Ppm=[]
        Pmp=[]
        Pmm=[]
        for k in range(len(self.detuning_list)):
            Bell_correlator.append(((np.abs(A[k]))**2+(np.abs(D[k]))**2-(np.abs(B[k]))**2-(np.abs(C[k]))**2)/((np.abs(A[k]))**2+(np.abs(D[k]))**2+(np.abs(B[k]))**2+(np.abs(C[k]))**2))
            Pmm.append((np.abs(A[k]))**2)
            Ppp.append((np.abs(D[k]))**2)
            Pmp.append((np.abs(B[k]))**2)
            Ppm.append((np.abs(C[k]))**2)

        return( Bell_correlator, Pmm, Ppp, Pmp, Ppm )

    def plot_pulse(self):
        pulse_abs = []
        for k in range(len(self.time)):
            pulse_abs.append(np.abs(self.Bell_pulse(self.time[k])))
        plt.figure(figsize=(10,7))
        plt.xlabel("Time (ms)")
        plt.ylabel(r"$|\Omega_R|$ (kHz)")
        plt.grid(True)
        plt.plot(self.time * 1e3, np.array(pulse_abs)/(2*pi)*1e-3)
        plt.show()


    # Detuning plot
    def plot_boxes(self,plt):
        if self.show_boxes is True:
            for k in range(self.boxes_number+1):
                plt.axvline(x=-((2*k+1)/2)*self.box_width/(2*pi) *1e-3,color='black',linestyle='--',linewidth=1)
                plt.axvline(x=+((2*k+1)/2)*self.box_width/(2*pi) *1e-3,color='black',linestyle='--',linewidth=1)
            plt.axvspan(-self.box_width/(2*pi) *(1e-3)/2-self.boxes_number*self.box_width/(2*pi) *(1e-3), -self.box_width/(2*pi) *(1e-3)/2, alpha=0.2, color='red')
            plt.axvspan(self.box_width/(2*pi) *(1e-3)/2, self.box_width/(2*pi) *(1e-3)/2+self.boxes_number*self.box_width/(2*pi) *(1e-3),alpha=0.2, color='green')


    def plot_probabilities(self,results):
        plt.figure(figsize=(10,7))
        plt.grid(True)
        self.plot_boxes(plt)
        plt.xlabel("Detuning (kHz)")
        plt.ylabel("Probability")
        ( Bell_correlator, Pmm, Ppp, Pmp, Ppm ) = results
        plt.plot(self.detuning_list / (2*pi) *1e-3, Pmm, label = r'$P_{--}$')
        plt.plot(self.detuning_list / (2*pi) *1e-3, Ppp, linestyle='--', label = r'$P_{++}$')
        plt.plot(self.detuning_list / (2*pi) *1e-3, Ppm, label = r'$P_{+-}$')
        plt.plot(self.detuning_list / (2*pi) *1e-3, Pmp, linestyle='--', label = r'$P_{-+}$')
        plt.title("Joint probabilities of detection")
        plt.legend()
        plt.show()

    def plot_correlator(self,results):
        plt.figure(figsize=(10,7))
        plt.grid(True)
        self.plot_boxes(plt)
        plt.xlabel("Detuning (kHz)")
        plt.ylabel("E")
        plt.xlim(0,np.max(self.detuning_list / (2*pi) *1e-3))
        for k in range(len(results)):
            ( Bell_correlator, Pmm, Ppp, Pmp, Ppm, label_curve ) = results[k]
            plt.plot(self.detuning_list / (2*pi) *1e-3, Bell_correlator,  label = label_curve)
        plt.title("Bell correlator")
        plt.legend()
        #plt.axhline(y=np.sqrt(2)/2,linestyle='--', color = 'black')
        plt.show()

# TOOLS
    def display(self):
        print(self.t_start_mirror)
        print(self.t_stop_mirror)
        print(self.splitter.pulse(0.0005))


