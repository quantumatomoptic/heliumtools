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

    def get_potential(self, X: torch.tensor, Y: torch.tensor, Z: torch.tensor, time: float = None):
        return torch.zeros_like(X)


class Trap(LinearPotential):
    """Harmonic trapping potential

    Args:
        omegax (Union[float, Callable]): The frequency along the x axis of the harmonic oscillator in Hz. It can be set to be either a constant or a function of time.
        omegay (Union[float, Callable]): The frequency along the y axis of the harmonic oscillator in Hz. It can be set to be either a constant or a function of time.
        omegaz (Union[float, Callable]): The frequency along the z axis of the harmonic oscillator in Hz. It can be set to be either a constant or a function of time.
    """

    def __init__(self, omegax: Union[float, Callable], omegay: Union[float, Callable], omegaz: Union[float, Callable]):
        super().__init__()

        self.omegax = omegax
        self.omegay = omegay
        self.omegaz = omegaz

    def on_propagation_begin(self):
        self.is_time_dependent = any_time_dependent_variable(self.omegax, self.omegay, self.omegaz)

        self._omegax = time_dependent_variable(self.omegax)
        self._omegay = time_dependent_variable(self.omegay)
        self._omegaz = time_dependent_variable(self.omegaz)

    def get_potential(self, X: torch.tensor, Y: torch.tensor, Z: torch.tensor, time: float = None):
        return 0.5*(1/self.gas.adim_pulse)**2*((self._omegax(time)*X)**2+(self._omegay(time)*Y)**2+(self._omegaz(time)*Z)**2)


# --- Non linear potentials ---
class Contact(NonLinearPotential):
    """Contact interactions potential

    Args:
        a_s (float): The scattering length in units of the Bohr radius.
        a_orth (float): The normalization parameter for the axis.
    """

    def __init__(self, a_s: float = 100, a_orth: float = 1e-6):
        super().__init__()

        self.a_s = a_s
        self.a_orth = a_orth

    def on_propagation_begin(self):
        self._a_s = self.a_s*spconsts.codata.value("Bohr radius")
        self._g = 4*np.pi*self.gas.N_particles*self._a_s/self.a_orth

    def potential_function(self, X: torch.tensor, Y: torch.tensor,  Z: torch.tensor, psi: torch.tensor, time: float = None):
        return self._g*torch.abs(psi)**2