from __future__ import annotations
from typing import Union, Callable

import scipy.constants as spconsts
import numpy as np
import torch

from ..utils.potentials import LinearPotential, NonLinearPotential, any_time_dependent_variable, time_dependent_variable


# --- Linear potentials ---
class Zero(LinearPotential):
    """Zero potential. It is equivalent to not applying any potential at all.
    """

    def __init__(self):
        super().__init__()

    def get_potential(self, X: torch.tensor, time: float = None):
        return torch.zeros_like(X)


class Trap(LinearPotential):
    """Harmonic trapping potential

    Args:
        omegax (Union[float, Callable]): The frequency along the x axis of the harmonic oscillator in Hz. It can be set to be either a constant or a function of time.
        omegay (Union[float, Callable]): The frequency along the y axis of the harmonic oscillator in Hz. It can be set to be either a constant or a function of time.
    """

    def __init__(self, omegax: Union[float, Callable]):
        super().__init__()

        self.omegax = omegax

    def on_propagation_begin(self):
        self.is_time_dependent = any_time_dependent_variable(self.omegax)

        self._omegax = time_dependent_variable(self.omegax)

    def get_potential(self, X: torch.tensor, time: float = None):
        return 0.5*(1/self.gas.adim_pulse)**2*(self._omegax(time)*X)**2

class GaussianTrap(LinearPotential):
    """Harmonic trapping potential

    Args:
        omegax (Union[float, Callable]): The frequency along the x axis of the harmonic oscillator in Hz. It can be set to be either a constant or a function of time.
        omegay (Union[float, Callable]): The frequency along the y axis of the harmonic oscillator in Hz. It can be set to be either a constant or a function of time.
    """

    def __init__(self, omegax: Union[float, Callable]):
        super().__init__()

        self.omegax = omegax

    def on_propagation_begin(self):
        self.is_time_dependent = any_time_dependent_variable(self.omegax)

        self._omegax = time_dependent_variable(self.omegax)

    def get_potential(self, X: torch.tensor, time: float = None):
        return 0.5*(1/self.gas.adim_pulse)**2*(self._omegax(time)*X)**2


class Lattice(LinearPotential):
    """Lattice potential

    Args:
        kLatt (float) : The lattice wavevector.
        V0 (Union[float, Callable]): The lattice depth in units of the recoil energy. It can be set to be either a constant or a function of time.
        theta (float): The angle of the lattice in the 2D plane.
        phi (Union[float, Callable]): The phase of the lattice. Can be time dependent, e.g, for Bragg diffraction
    """

    def __init__(self, kLatt: float = 0.0, V0: Union[float, Callable] = 0, phi: Union[float, Callable] = 0):
        super().__init__()

        self.kLatt = kLatt
        self.V0 = V0
        self.phi = phi

    def on_propagation_begin(self):
        # compute recoil energy
        self.Er = 0.5 * (spconsts.hbar*self.kLatt)**2 / self.gas.mass
        # prepare V0 and phi
        self.is_time_dependent = any_time_dependent_variable(self.V0, self.phi)
        self._V0 = time_dependent_variable(self.V0)
        self._phi = time_dependent_variable(self.phi)

    def get_potential(self, X: torch.tensor, time: float = None):
        phase = self.kLatt*X*self.gas.adim_length - self._phi(time)
        arg = 1 + torch.cos(phase)
        return self._V0(time)*self.Er/(spconsts.hbar*self.gas.adim_pulse)*arg


class SquareBox(LinearPotential):
    """Square box potential

    Args:
        V (float): The depth of the box.
        D (float): The size of the box.
    """

    def __init__(self, V: float, D: float):
        super().__init__()
        self.V = V
        self.D = D

    def on_propagation_begin(self):
        self._V = self.V/(spconsts.hbar*self.gas.adim_pulse)
        self._D = self.D/self.gas.adim_length
        self._box = self._V * (1-torch.heaviside(-torch.abs(self.gas.X)+self._D/2, torch.ones_like(self.gas.X)))

    def get_potential(self, X: torch.tensor, time: float = None):
        return self._box


# --- Non linear potentials ---
class Contact(NonLinearPotential):
    """Contact interactions potential

    Args:
        a_s (float): The scattering length in units of the Bohr radius.
        a_orth (float): The renormalization parameter for the scattering length to account for the missing third dimension.
    """

    def __init__(self, a_s: float = 100, a_orth: float = 1e-6):
        super().__init__()

        self.a_s = a_s
        self.a_orth = a_orth

    def on_propagation_begin(self):
        self._a_s = self.a_s*spconsts.codata.value("Bohr radius")
        self._g = np.sqrt(8*np.pi)*self.gas.N_particles*self._a_s/self.a_orth

    def potential_function(self, X: torch.tensor, psi: torch.tensor, time: float = None):
        return self._g*torch.abs(psi)**2