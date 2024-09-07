from scipy.constants import k as kb
from scipy.constants import hbar, m_p, h
from heliumtools.misc.logger import getLogger, logging
import scipy
import numpy as np

log = getLogger(__name__)
# set the desired level of warning : DEBUG / INFO / WARNING
log.setLevel(logging.DEBUG)


class BEC:
    _asn1 = 1
    _a_s = 7e-9
    _m = 4 * m_p
    _g = 4 * np.pi * hbar**2 * _a_s / _m
    _omega_perp = 2 * np.pi * 1000  # Hz
    _omega_z = 2 * np.pi * 30  # Hz
    _temperature = 30 * 1e-9
    _c_s = None

    def __init__(self, **kwargs) -> None:
        self.__dict__.update(**kwargs)
        self.update_trap_parameters()

    ###########################
    ## DEFINITION OF VARIABLES
    ###########################
    @property
    def omega_perp(self):
        return self._omega_perp

    @omega_perp.setter
    def omega_perp(self, value):
        self._omega_perp = value
        self.update_trap_parameters()
        self.set_up_BEC_properties()

    @property
    def omega_z(self):
        return self._omega_z

    @omega_z.setter
    def omega_z(self, value):
        self._omega_z = value
        self.update_trap_parameters()
        self.set_up_BEC_properties()

    @property
    def asn1(self):
        return self._asn1

    @asn1.setter
    def asn1(self, value):
        self._asn1 = value
        self.set_up_BEC_properties()

    @property
    def temperature(self):
        return self._temperature

    @temperature.setter
    def temperature(self, value: float):
        """set the temperature of the gas, in K

        Parameters
        ----------
        value : float
            temperature of the gas, in K
        """
        self._temperature = value
        self.set_up_BEC_properties()

    ###########################
    ## DEFINITION OF METHODS
    ###########################
    def update_trap_parameters(self):
        """function that updates the trap parameters using omega_perp and omega_z values"""
        self._a_perp = np.sqrt(hbar / (self._m * self._omega_perp))
        self._a_z = np.sqrt(hbar / (self._m * self._omega_z))
        self._lamb = self._omega_z / self._omega_perp
        if self._lamb >= 1:
            msg = "Your BEC is a pancake and not a cigar. The transverse freq is "
            msg += f"{self.omega_perp/(2*np.pi):.3e} Hz and the longitudinal is {self.omega_z/(2*np.pi):.3e} Hz."
            log.warning(msg)

    def set_longitudinal_freq(self, freq):
        """Set the vertical (longitudinal) frequency of the trap of the BEC

        Parameters
        ----------
        freq : float
            longitudinal frequency of the trap in Hz (do not multiply by 2pi)
        """
        self.omega_z = 2 * np.pi * freq
        # self.update_trap_parameters()

    def set_transverse_freq(self, freq):
        """Set the horizontal frequency of the trap of the BEC

        Parameters
        ----------
        freq : float
            transverse frequency of the trap in Hz  (do not multiply by 2pi)
        """
        self.omega_perp = 2 * np.pi * freq
        # self.update_trap_parameters()

    def set_up_BEC_properties(self, asn1=10):
        """Update all BEC properties based on the value of asn1 at the center of the trap."""
        log.warning(
            "The BEC class is generic and does not contain specific model. Please instanciate a model."
        )

    ####################################
    ## DEFINITION OF ADVANCED METHODS ##
    ####################################
    def get_phonon_energy(self, v_ph: float):
        """Return the phonon energy of the Bogoliubov spectrum.

        v_ph : float
            speed at which one wants the thermal population.
        """
        k = self._m * v_ph / hbar
        epsilon_k = hbar**2 * k**2 / (2 * self._m)
        energy = np.sqrt(2 * self._mc2 * epsilon_k + epsilon_k**2)
        return energy

    def get_phonon_population(self, v_ph: float):
        """return the phonon population at speed v_ph of a homogeneous BEC

        Parameters
        ----------
        v_ph : float
            speed at which one wants the thermal population.
        """
        energy = self.get_phonon_energy(v_ph)
        n_th = 1 / (np.exp(energy / self._temperature / kb) - 1)
        return n_th

    def get_phonon_population_from_energy(self, energy: float):
        """return the phonon population of the BEC

        Parameters
        ----------
        energy : float
            phonon energy in SI.
        """
        n_th = 1 / (np.exp(energy / self._temperature / kb) - 1)
        return n_th

    def get_phonon_group_speed(self, v_ph: float):
        """return the phonon group speed velocity in the homogenous BEC limit

        Parameters
        ----------
        v_ph : float
            phonon speed in m/s (SI units).

        Returns
        -------
        float
            phonon group velocity in m/s (SI units)

        """
        k = self._m * v_ph / hbar
        epsilon_k = hbar**2 * k**2 / (2 * self._m)
        energy = np.sqrt(2 * self._mc2 * epsilon_k + epsilon_k**2)
        deriv_energysquare = 4 * self._mc2 * epsilon_k / k + 4 * epsilon_k**2 / k
        return 0.5 * deriv_energysquare / energy / hbar



    
class CigarShapeBEC(BEC):
    _c_s = None
        
    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self.set_up_BEC_properties() 



    def set_up_BEC_properties(self):
        """set up the BEC properties based on the value of asn1, in the Gaussian ansatz approximation.
        References:
        * Gerbier EPL (2004)
        * Petrov et al PRL (2000)
        * Petrov et al  J. Phys (2004)
        * Victor's thesis for notations.
        """
        self._mu0 = hbar * self._omega_perp * 2 * np.sqrt(self._asn1)
        self._radius_perp = 2 * self._a_perp * self._asn1**0.25
        # self._chi = np.sqrt(self._alpha**3 * (5 + self._alpha**2))
        self._c_s = np.sqrt(
            self._mu / 2
            / self._m
        )
        self._Nat = int(
            self._a_z**2
            / self._a_s
            / self._a_perp
            / 15
            * (self._alpha ** (1.5) * (self._alpha + 5))
        )  # atom number
        self._length = (
            np.sqrt(self._alpha * self._a_z**4 / self._a_perp**2) * 2
        )  # total lenght of bec
        self._n1p = self._asn1 / self._a_s  # max density

        self._radius = (
            self._a_perp * (1 + 4 * self._asn1) ** 0.25
        )  # bec transverse sigma
        self._mc2 = self._m * self._c_s**2  # mc² phonon energy
        self._xi = hbar / (self._m * self._c_s * np.sqrt(2))  # healing length
        self._g1D = self._g / (2 * np.pi * self._sigma0**2)
        self._gamma = self._m * self._g1D / hbar**2 / self._n1p  #
        # self._radial_0
        self._l_coh = hbar / np.sqrt(
            self._m * self._n1p * self._g1D
        )  # coherence lenght, Petrov and
        self._T_phi = (
            15
            * (hbar * self._omega_z) ** 2
            / (32 * (self._mu0 - self._omega_perp * hbar))
            * self._Nat
            / kb
        )
        self._length_phase = 2 * self._n1p * hbar**2 / kb / self._temperature / self._m
        # round(hbar/(m_he * (1e-3*cs) * np.sqrt(2))*1e6,1)
        # log.info("Instanciation faite.")




class Gaussian_BEC(BEC):
    _c_s = None

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self.set_up_BEC_properties()

    @property
    def c_s(self):
        return self._c_s

    @c_s.setter
    def c_s(self, val):
        self.set_sound_speed(val)

    def evaluate_1D_density(self, z=np.array([])) -> np.array:
        """return the BEC 1D dentisy evaluated at position z.

        Parameters
        ----------
        z : np.array, optional
            the point where the density should be evalutated, by default np.array([])

        Returns
        -------
        np.array
            the density in unit of asn1
        """
        radius = self._length / 2
        if len(z) == 0:
            z = np.linspace(-radius, radius, 200)

        return np.where(
            np.abs(z) < radius,
            self._alpha
            / 16
            * (1 - z**2 / radius**2)
            * (self._alpha * (1 - z**2 / radius**2) + 4),
            z * 0,
        )

    def set_sound_speed(self, cs):
        """Set the BEC parameter from the measured sound speed

        Parameters
        ----------
        cs : float
            BEC speed of sound (m/s, SI units)
        """
        mc2 = self._m * cs**2
        hbaromega = hbar * self._omega_perp
        asn1 = 0.5 * (mc2 / hbaromega) ** 2 * (1 + np.sqrt(1 + (hbaromega / mc2) ** 2))
        self.asn1 = asn1

    def set_sound_speed_from_parametric_resonance(self, v_ph):
        """From the Bogoliubov phonon speed of sound, it compute the speed of sound from which
        the class infers the asn1 parameter and set the BEC properties.

        Parameters
        ----------
        v_ph : float
            The phonon speed of the resonance Bogoliubov pairs (m/s, SI units)
        """
        k = self._m / hbar * v_ph
        cs = float(
            np.sqrt(self._omega_perp**2 - (hbar * k**2 / (2 * self._m)) ** 2) / k
        )  # m/s
        self.set_sound_speed(cs)

    def evaluate_local_speed_of_sound(self, z=np.array([])) -> np.array:
        """return the BEC 1D speed of sound in the LDA evaluated at position z.

        Parameters
        ----------
        z : np.array, optional
            the point where the cs should be evaluated, by default np.array([])

        Returns
        -------
        np.array
            the speed of sound in the LDA
        """
        asn1_z = self.evaluate_1D_density(z)
        return (
            self._c_s
            * np.sqrt(asn1_z / np.sqrt(1 + 4 * asn1_z))
            / (np.sqrt(self._asn1 / np.sqrt(1 + 4 * self._asn1)))
        )

    def get_bogoliubov_spectrum(self, nu=np.arange(1, 10)) -> np.array:
        """return excitation spectrum of longitudinal excitations.
        See Schemmer's thesis page 27 or Pitaevskii & Stringari

        Parameters
        ----------
        nu : _type_, optional
            _description_, by default np.arange(1, 10 )

        Returns
        -------
        np.array
            _description_
        """

        return np.sqrt(nu * (nu + 1) / 2) * self._omega_z

    def get_homogeneous_bogoliubov_spectrum(self, nu=np.arange(1, 10)):
        """return the excitation spectrum for all  k = 2 nu / L, where L is the BEC length and nu is the excitation number."""
        k = nu / (self._length / 2)
        epsilon_k = hbar**2 * k**2 / (2 * self._m)
        omega = np.sqrt(self._mc2 * epsilon_k + epsilon_k**2) / hbar
        return omega

    def get_excitation_form(
        self,
        z=np.array([]),
        nu=1,
    ) -> np.array:
        """return the density Bogoliubov excitation wavefunction of the excitation number nu. See Schemmer et al. PRA98 (2018) for formulae. See also Petrov PRL (2000).


        Parameters
        ----------
        z : in situ position, optional
            relative position of the trap, by default np.array([])
        nu : int, optional
            number of the excitation, by default 1

        Returns
        -------
        np.array
            density fluctuations as a function of z.
        """
        if len(z) == 0:
            z = np.linspace(-1, 1, 200)
        n_p = self._asn1 / self._a_s  # peak density
        if type(nu) != int:
            log.error(
                f"GaussianBEC: excitation number should be an integer, it is {type(nu)}. Return BEC mode."
            )
            return self.evaluate_1D_density(z)
        legendre_poly = scipy.special.eval_legendre(nu, z)
        density_fluc = (
            np.sqrt(2 * nu + 1)
            / self._length
            * (nu * (nu + 1)) ** 0.25
            * (hbar**2 * n_p / self._m / self._g) ** 0.25
            * legendre_poly
            * self._a_s
        )
        return density_fluc

    def set_up_BEC_properties(self):
        """set up the BEC properties based on the value of asn1, in the Gaussian ansatz approximation.
        References:
        * Gerbier EPL (2004)
        * Petrov et al PRL (2000)
        * Petrov et al  J. Phys (2004)
        * Victor's thesis for notations.
        """
        self._mu0 = hbar * self._omega_perp * np.sqrt(1 + 4 * self._asn1)
        self._sigma0 = self._a_perp  * (1 + 4 * self._asn1)**(0.25)
        self._alpha = 2 * (np.sqrt(1 + 4 * self._asn1) - 1)
        self._chi = np.sqrt(self._alpha**3 * (5 + self._alpha**2))
        self._c_s = np.sqrt(
            self.asn1
            / self._m
            * 2
            * hbar
            * self._omega_perp
            / np.sqrt(1 + 4 * self._asn1)
        )
        self._Nat = int(
            self._a_z**2
            / self._a_s
            / self._a_perp
            / 15
            * (self._alpha ** (1.5) * (self._alpha + 5))
        )  # atom number
        self._length = (
            np.sqrt(self._alpha * self._a_z**4 / self._a_perp**2) * 2
        )  # total lenght of bec
        self._n1p = self._asn1 / self._a_s  # max density

        self._radius = (
            self._a_perp * (1 + 4 * self._asn1) ** 0.25
        )  # bec transverse sigma
        self._mc2 = self._m * self._c_s**2  # mc² phonon energy
        self._xi = hbar / (self._m * self._c_s * np.sqrt(2))  # healing length
        self._g1D = self._g / (2 * np.pi * self._sigma0**2)
        self._gamma = self._m * self._g1D / hbar**2 / self._n1p  #
        # self._radial_0
        self._l_coh = hbar / np.sqrt(
            self._m * self._n1p * self._g1D
        )  # coherence lenght, Petrov and
        self._T_phi = (
            15
            * (hbar * self._omega_z) ** 2
            / (32 * (self._mu0 - self._omega_perp * hbar))
            * self._Nat
            / kb
        )
        self._length_phase = 2 * self._n1p * hbar**2 / kb / self._temperature / self._m
        # round(hbar/(m_he * (1e-3*cs) * np.sqrt(2))*1e6,1)
        # log.info("Instanciation faite.")

    def show_bec_info(self):
        msg = "Gas parameters in Gaussian ansatz \n" + "=" * 30 + "\n"
        msg += (
            "1D parameters is asn1 = {:.1e} (asN/L={:.1e}) and chi = {:.1e}\n".format(
                self._asn1, self._Nat / self._length * self._a_s, self._chi
            )
        )
        msg += "It is a trap of {:.0f} bosons with frequency ".format(self._Nat)
        msg += "{:.1f} Hz and\n{:.2f} kHz.".format(
            self._omega_z / 2 / np.pi, self._omega_perp / 2 / np.pi * 1e-3
        )
        msg += "The peak atomic density is {:.2e} at/cm^3\n".format(
            self._n1p / self._radius**2 / 1e6
        )
        msg += "The BEC length is {:.0f} µm with a radius of {:.1f} µm.\n".format(
            self._length * 1e6, self._radius * 1e6
        )
        msg += "It can be compared to the correlation radius {:.1f} µm \n".format(
            self._length_phase * 1e6
        )
        msg += "The interparticle mean distance is {:.1f} nm\n".format(
            (self._length * self._sigma0**2 * 4 / self._Nat) ** (1 / 3) * 1e9
        )
        msg += "that should be compare to 7 nm scattering length.\n"
        msg += (
            "The healing length is {:.3f} µm and the sound speed {:.3f} mm/s.".format(
                self._xi * 1e6, self._c_s * 1e3
            )
        )
        msg += "\nThe lieb liniger parameter is gamma = {:.1e}".format(self._gamma)
        msg += "\nand the ratio kbT/mc2 = {:.2f}.".format(
            kb * self._temperature / self._mc2
        )
        # msg += "\nTphi={:.2f} nK".format(self._T_phi*1e9)
        log.info(msg)


if __name__ == "__main__":
    bec = Gaussian_BEC()
    bec.set_transverse_freq(1.8e3)
    phonon_speed = 17/2/1000
    bec.set_longitudinal_freq(30)
    bec.set_sound_speed_from_parametric_resonance(phonon_speed)
    bec.temperature = 40*1e-9

    bec.show_bec_info()
    bec.get_phonon_population(phonon_speed)


