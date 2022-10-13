#!/usr/bin/env python
# -*- mode:Python; coding: utf-8 -*-

# ----------------------------------
# Created on the Tue Oct 11 2022 by Victor, based on the code of Charlie
#
# Developped by Charlie, Victor
#
# Last (big) change on the ... by ...
#
# Copyright (c) 2022 - Helium1@LCF
# ----------------------------------
#
"""
Content of lattice.py

Orientation des axes : la gravité est orientée vers le bas : le faisceau up monte et le faisceau down descend. 





"""
import scipy.constants as cst
from heliumtools.units import u, const
from heliumtools.atom import Heliumqunits
import numpy as np
import warnings
from math import pi
from numpy import linalg as LA
from heliumtools import qunits
from scipy import interpolate


class Lattice:
    """
    Lattice Attributs

    up_frequency ; down_frequency
    up_power ; down_power
    up_waist ; down_waist
    theta : angle between the two beams (radians)
    wavelength : lattice wavelength
    atom : class atom (helium by default)
    up_intensity ; down_intensity : intensité du réseau
    detuning : detuning up - down du réseau (qunits)
    k : vecteur d'onde du réseau (qunits)
    v : vitesse du réseau (qunits)
    E : énergie du réseau h^2 k^2 / 2m
    V0 : profondeur du réseau
    """

    def __init__(self, **kwargs):
        self.up_frequency = 200 * u.MHz
        self.down_frequency = 200.100 * u.MHz
        self.up_power = 85 * u.mW
        self.down_power = 85 * u.mW
        self.up_waist = 200 * u.um
        self.down_waist = 200 * u.um
        self.theta = 166 / 360 * 2 * pi
        self.wavelength = 1064 * u.nm
        self.atom = Heliumqunits()
        self.Cjmatrix_size = 41  # see Dussarat PhD thesis page 45
        self.__dict__.update(kwargs)
        self.build_lattice_properties()

    def build_lattice_properties(self):
        """Construit les propriétés de la classe Lattice à partir des arguments initialisés. Cette fonction est appelée en fin d'initialisation."""
        self.up_intensity = 2 * self.up_power / self.up_waist**2 / pi
        self.down_intensity = 2 * self.down_power / self.down_waist**2 / pi
        self.detuning = 2 * pi * (self.up_frequency - self.down_frequency)
        self.k = np.sin(self.theta / 2) * 2 * pi / self.wavelength
        self.v = self.detuning / 2 / self.k
        self.E = const.hbar**2 * self.k**2 / 2 / self.atom.mass
        if self.down_intensity != self.up_intensity:
            warnings.warn(
                "Up and down arms do not have the same intensity. I will take the mean to define the intensity. "
            )
        self.intensity = np.sqrt(self.down_intensity * self.up_intensity)
        self.V0 = const.hbar * self.atom.atomic_width**2 * self.intensity
        ## TO BE CHECKED BY SOMEONE ELSE
        alpha = self.atom.get_alpha(self.wavelength)  # Polarisability
        self.V0 = np.abs(0.5 / const.eps0 / const.c * alpha * self.intensity)

    def check_attributs(self):
        """fonction appelée après l'itnitialisation pour vérifier tous les attributs."""
        self.check_Cjmatrixsize()

    def set_Cjmatrixsize(self, size):
        self.Cjmatrix_size = int(size)
        self.check_Cjmatrixsize()

    def check_Cjmatrixsize(self):
        self.Cjmatrix_size = int(self.Cjmatrix_size)
        if self.Cjmatrix_size < 3:
            self.Cjmatrixsize == 41
            warnings.warn(
                f"pThe Cj matrix size must be big enough. Please check Dussarat PhD thesis page 45."
            )
        if self.Cjmatrix_size % 2 == 0:
            warnings.warn(
                f"The Cj matrix size must be odd. It was set to {self.Cjmatrixsize+1}"
            )
            self.Cjmatrixsize += 1

    def get_bands_structure(self, n=0, precision=1001, max=1):
        """Retourne l'énergie de la bande numéro n.

        Parameters
        ----------
        n : int, optional
            numéro max de la bande, 0 bande fondamentale, 1 prmièer bande excitée etc..., by default 0
        precision : int, optional
            nombre de quasi momentum, par défaut 1001 --> doit être impair

        Returns
        -------
        2 lists of length precision with quasimomentum and corresponding energy.
        """
        if precision < 0:
            warnings.warn(f"precision was negative. get_band() returns empty lists.")
            return ([], [])
        if precision % 2 == 0:
            precision += 1
            warnings.warn(
                f"precision was even and must be odd : precision set to {precision}"
            )
        # On définit ensuite une taille de matrice pas trop grosse par rapport à la précision voulue.
        if n < 3:
            self.set_Cjmatrixsize(5)
        else:
            self.set_Cjmatrixsize(n)
        q_list = np.linspace(-1, 1, precision) * self.k
        energy_list = np.zeros((self.Cjmatrix_size, precision))
        for i, q in enumerate(q_list):
            matrix = self._generate_matrix(q)
            w, v = LA.eig(matrix)
            w.sort()
            energy_list[:, i] = w
        energy_list *= self.E  # convert energy to its real unit.
        return (q_list, energy_list)

    ######################################################
    ### INTERNAL FUNCTIONS FOR CALCULATIONS
    ######################################################

    def _generate_matrix(self, q):
        """Génère la mtrice des coefficients Cj permettant d'avoir fonctions d'ondes et énergies dans le réseau.

        Parameters
        ----------
        q : qunits inverse of length
            vecteur d'onde pour lequel on veut déterminé la matrice.

        Returns
        -------
        _type_
            _description_
        """
        self.check_Cjmatrixsize()
        matrix = np.zeros((self.Cjmatrix_size, self.Cjmatrix_size))
        nbsim = int((self.Cjmatrix_size - 1) / 2)
        for k in range(self.Cjmatrix_size):
            j = k - nbsim
            matrix[k, k] = (q / self.k + 2 * j) ** 2 + self.V0 / (2 * self.E)
            if k != self.Cjmatrix_size - 1:
                matrix[k, k + 1] = -self.V0 / (4 * self.E)
                matrix[k + 1, k] = -self.V0 / (4 * self.E)
        return matrix

    def momentum_to_quasimomentum(self, k, units=False):
        if units:
            adim_k = (k / self.k).to(u.dimensionless)
            adim_q = self.momentum_to_quasimomentum(adim_k, units=False)
            return adim_q * self.k
        else:
            return k - np.sign(k) * 2 * ((k + 1) // 2)

    def get_pairs_quasi_momentum(self, bec_speed=0 * u.mm / u.s, precision=201):
        bec_quasi_momentum = (
            self.atom.speed_to_momentum(bec_speed - self.v) / self.k
        ).to(u.dimensionless)
        ## D'abord, il faut trouver la relation de disperion de la bande fondamentale
        n_bands = 5
        quasimomentum, energies = self.get_bands_structure(
            n=n_bands, precision=precision
        )
        # On récupère ensuite seulement la bande fondamentale. On traite maintenant le problème SANS DIMENSION car concaténer des arrays avec des dimensions pose problème.
        fondam = (energies[0, :] / self.E).to(u.dimensionless)
        fondam = np.concatenate((fondam, fondam, fondam))
        momentum = (quasimomentum / self.k).to(u.dimensionless)
        momentum = np.concatenate((momentum - 2, momentum, momentum + 2))
        ### RECHERCHE DE LA PAIRE
        ## On définit une fonction E (l'énergie du fondamental) qui interpole notre résultat
        E = interpolate.interp1d(momentum, fondam)
        # On recherche maintenant les paires (recherche dummy de racine)
        def energy_conservation(q):
            return 2 * E(bec_quasi_momentum) - E(q) - E(2 * bec_quasi_momentum - q)

        # Now we want to fin the zero of the function (I cannot use fsolve because the function is not defined evrywhere and it brake the programm.)
        q = bec_quasi_momentum - 1
        step = 1 / 10000.0
        if precision > 10000:
            step = 1 / precision
        timeout = False
        while (
            q < bec_quasi_momentum + 1
            and energy_conservation(q) * energy_conservation(q + step) > 0
        ):
            q += step
        if timeout:
            q1 = bec_quasi_momentum
            q2 = bec_quasi_momentum
        else:
            q1 = q
            q2 = 2 * bec_quasi_momentum - q

        ## CONVERSION MOMENTUM --> QUASIMOMENTUM

        return (self.momentum_to_quasimomentum(q1), self.momentum_to_quasimomentum(q2))
