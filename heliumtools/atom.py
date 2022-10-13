# -*- coding: utf-8 -*-
"""
Author   : alex
Created  : 2020-10-13 17:06:38
Modified : 2020-10-19 17:19:33

Comments : contains atomic data (for Helium)
"""
# -- imports
import scipy.constants as csts
from heliumtools.units import u, const

import numpy as np


class Heliumqunits:
    """
    Atomic data & calculations for Helium using qunits.
    """

    def __init__(self, **kwargs):
        """
        Object initialization, sets parameters
        """
        # -- initialize default settings
        # physical parameters
        self.mass = 4.002602 * const.u  # u.kg
        self.atomic_wavelength = 1083 * u.nm
        self.atomic_width = 2 * np.pi * 1.62 * u.MHz
        self.atomic_omega = 2 * np.pi * const.c / self.atomic_wavelength
        self.scattering_length = 7.512 * u.nm  # m, s-wave scattering length
        # magnetic
        self.lande_g_factor = 2
        self.zeeman_state = -1

    def get_alpha(self, wavelength=1550 * u.nm, unit="SI"):
        """
        computes the polarizability.
        input  : wavelength = laser wavelength (qunits)
        output : polarizability (qunits)
        """
        # parameters
        omega_L = 2 * np.pi * const.c / wavelength
        omega_0 = self.atomic_omega
        Gamma = self.atomic_width

        # compute
        Delta = 1 / (1 / (omega_0 - omega_L) + 1 / (omega_0 + omega_L))
        alpha = 3 * np.pi * const.eps0 * const.c**3 * Gamma / omega_0**3 / Delta
        return alpha

    def get_scattering_rate(self, intensity, wavelength=1550 * u.nm):
        """
        computes the scattering rate at given intensity & detuning
        input  : intensity = laser intensity (in qunits)
                 wavelength = laser wavelength (in qunits)
        output : scattering rate (in qunits)
        """
        # parameters
        omega_L = 2 * np.pi * const.c / wavelength
        omega_0 = self.atomic_omega
        Gamma = self.atomic_width

        # compute
        Delta = 1 / (1 / (omega_0 - omega_L) + 1 / (omega_0 + omega_L))
        Gamma_sc = (
            3
            * np.pi
            * const.c**2
            / 2
            / const.hbar
            / omega_0**3
            * (omega_L / omega_0) ** 3
            * (Gamma / Delta) ** 2
            * intensity
        )
        return Gamma_sc

    def speed_to_momentum(self, v):
        return self.mass * v / const.hbar

    def convert_speed_to_lattice_momentum(
        self, v, wavelength=1064 * u.nm, theta=166 * 2 * np.pi / 360
    ):
        """Convert atom speed (mm/s) to momentum in hbar unit (SI)

        Parameters
        ----------
        v : float
            speed of the atom in qunits
        wavelength : float, optional
            wavelenght of the lattice (in nanometers), by default 1064
        theta : int, optional
            angle in degrees beween the two beams, by default 166

        Returns
        -------
        relative :float
            speed compare to lattice momentum
        """
        klatt = 2 * np.pi / wavelength * np.sin(theta)

        return self.mass * v / (const.hbar * klatt)


# -- define atom data dictionnary
class Helium:
    """
    Atomic data & calculations for Helium
    """

    def __init__(self, **kwargs):
        """
        Object initialization, sets parameters
        """
        # -- initialize default settings
        # physical parameters
        self.mass = 4.002602 * csts.m_u  # kg
        self.atomic_wavelength = 1083e-9  # m
        self.atomic_width = 2 * np.pi * 1.62e6  # s^-1
        self.atomic_omega = 2 * np.pi * csts.c / self.atomic_wavelength
        self.scattering_length = 7.512e-9  # m, s-wave scattering length
        # magnetic
        self.lande_g_factor = 2
        self.zeeman_state = -1

    def get_alpha(self, wavelength=1550e-9, unit="SI"):
        """
        computes the polarizability.
        input  : wavelength = laser wavelength (in meters)
                 unit = 'SI' (default) or 'au' (atomic units)
        output : polarizability (in SI / au)
        """
        # parameters
        omega_L = 2 * np.pi * csts.c / wavelength
        omega_0 = self.atomic_omega
        Gamma = self.atomic_width

        # compute
        Delta = 1 / (1 / (omega_0 - omega_L) + 1 / (omega_0 + omega_L))
        alpha = 3 * np.pi * csts.epsilon_0 * csts.c**3 * Gamma / omega_0**3 / Delta
        if unit == "SI":
            return alpha
        else:
            return alpha / 1.649e-41

    def get_scattering_rate(self, intensity, wavelength=1550e-9):
        """
        computes the scattering rate at given intensity & detuning
        input  : intensity = laser intensity (in W/m^2)
                 wavelength = laser wavelength (in meters)
        output : scattering rate (in s^-1)
        """
        # parameters
        omega_L = 2 * np.pi * csts.c / wavelength
        omega_0 = self.atomic_omega
        Gamma = self.atomic_width

        # compute
        Delta = 1 / (1 / (omega_0 - omega_L) + 1 / (omega_0 + omega_L))
        Gamma_sc = (
            3
            * np.pi
            * csts.c**2
            / 2
            / csts.hbar
            / omega_0**3
            * (omega_L / omega_0) ** 3
            * (Gamma / Delta) ** 2
            * intensity
        )
        return Gamma_sc

    def convert_speed_to_lattice_momentum(self, v, wavelength=1064.0, theta=166):
        """Convert atom speed (mm/s) to momentum in hbar unit (SI)

        Parameters
        ----------
        v : float
            speed of the atom in millimeter per second
        wavelength : float, optional
            wavelenght of the lattice (in nanometers), by default 1064
        theta : int, optional
            angle in degrees beween the two beams, by default 166

        Returns
        -------
        relative :float
            speed comapre to lattice momentum
        """
        klatt = 2 * np.pi / (wavelength * 1e-9) * np.sin(2 * np.pi * theta / 360)

        return 0.001 * self.mass * v / (csts.hbar * klatt)


# -- TESTS
if __name__ == "__main__":
    from laser import intensity_gauss

    # init object
    he = Helium()
    # polarizability
    print(he.get_alpha(unit="SI"))
    print(he.get_alpha(unit="au"))
    # scattering rate
    P = 5.6 * 1.2  # Watts
    w = 135e-6  # meters
    I0 = intensity_gauss(0, 0, w, P)
    print(he.get_scattering_rate(I0) * u.V / u.cm)
