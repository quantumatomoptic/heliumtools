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
from scipy.optimize import fsolve


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
        self.down_frequency = (
            200.120 * u.MHz
        )  # 90 kHz de différence donne 0.5klatt et 180 kHz donne 1klat
        self.up_power = 85 * u.mW
        self.down_power = 85 * u.mW
        self.up_waist = 585 * u.um
        self.down_waist = 585 * u.um
        self.theta = 166 / 360 * 2 * pi
        self.wavelength = 1064 * u.nm
        self.atomic_density = 1.3e13 / ((1 * u.cm) ** 3)
        self.atom = Heliumqunits()
        self.Cjmatrix_size = 41  # see Dussarat PhD thesis page 45
        self.clebsch = 1  # 1 / np.sqrt(3)
        self.__dict__.update(kwargs)
        self.build_lattice_properties()

    def upddate_property(self, **kwargs):
        """permet de mofifier des attributs de la classe lattice.  On update ensuite les propriétés qui découlent de la classe.
        Attention donc à ne pas modifier une propriété qui serait ensuite redéfinie dans la méthode build_lattice_properties (comme le detuning, issu des fréquences up et down du réseau).
        """
        self.__dict__.update(kwargs)
        self.build_lattice_properties()

    def build_lattice_properties(self):
        """Construit les propriétés de la classe Lattice à partir des arguments initialisés. Cette fonction est appelée en fin d'initialisation."""
        self.up_intensity = 2 * self.up_power / self.up_waist**2 / pi
        self.down_intensity = 2 * self.down_power / self.down_waist**2 / pi
        self.detuning = 2 * pi * (self.up_frequency - self.down_frequency)
        self.k = np.sin(self.theta / 2) * 2 * pi / self.wavelength
        self.v = self.detuning / 2 / self.k
        self.v_in_klatt = (self.atom.speed_to_momentum(self.v) / self.k).to(
            u.dimensionless
        )
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
        rabi_up = (
            self.atom.atomic_width
            * np.sqrt(self.up_intensity / self.atom.i_sat / 2)
            * self.clebsch
        )
        rabi_down = (
            self.atom.atomic_width
            * np.sqrt(self.down_intensity / self.atom.i_sat / 2)
            * self.clebsch
        )
        delta = (
            2 * pi * const.c * (1 / self.wavelength - 1 / self.atom.atomic_wavelength)
        )

        self.rabi_frequency = (
            rabi_up * rabi_down / delta / 2 / (2 * pi)
        )  # On exprime en Hz
        self.rabi_frequency2 = (
            2 * self.V0 / const.hbar / (2 * pi)
        )  # page 191 Quentin's thesis

    def show_lattice_properties(self, all=True):
        print("--------------------------")
        print(f"Lattice depth V0/E_latt = {self.V0/self.E}")
        print(f"Detuning δ/k_latt = {self.atom.speed_to_momentum(self.v) / self.k}")
        print(f"Rabi frequency = {self.rabi_frequency}")
        print("--------------------------")
        if all is True:
            for number, element in enumerate(self.__dict__):
                print(f"{element} : {self.__dict__[element]}")


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
            matrix = self._generate_Cjmatrix(q)
            w, v = LA.eig(matrix)
            w.sort()
            energy_list[:, i] = w
        energy_list *= self.E  # convert energy to its real unit.
        return (q_list, energy_list)

    def momentum_to_quasimomentum(self, k, units=False):
        """Méthode retournant la quasiimpulsion d'une particule étant donnée son impulsion. /!\ Si k est qunits d'inverse de longueur, il faut absolument le sépcifier en mettant unit à True.

        Parameters
        ----------
        k : sans dimension ou inverse de longueur (--> mais préciser dans ce cas qu'il a une unité en mettant units = True !!!)
            vecteur d'onde (array ou scalaire)

        units : bool, optional
            Si k a la dimension de l'inverse d'une longueur, mettre units =  True

        Returns
        -------
        quasi momentum dans la même unité que l'input.
        """
        if units:
            adim_k = (k / self.k).to(u.dimensionless)
            adim_q = self.momentum_to_quasimomentum(adim_k, units=False)
            return adim_q * self.k
        else:
            return k - 2 * ((k + 1) // 2)

    def get_pairs_quasi_momentum(self, bec_speed=0 * u.mm / u.s):
        # /!\/!\/!\
        # Dans la suite de la fonction, on travaille sans dimension !!
        # /!\/!\/!\
        bec_quasi_momentum = (
            self.atom.speed_to_momentum(bec_speed - self.v) / self.k
        ).to(u.dimensionless)
        if np.abs(bec_quasi_momentum) >= 1:
            warnings.warn(
                f"The BEC quasimomentum is greater than 1 ({bec_quasi_momentum}): this is not taken account in the model. Please modify me. "
            )
            return (bec_quasi_momentum, -1, 1)
        if np.abs(bec_quasi_momentum) < 0.5:
            warnings.warn(
                f"The BEC quasimomentum seems too low to creat pairs. Trying anyway but result might be false."
            )
        ## En faisant quelques test, je me rends compte que si q0 > 0, il fatu que la taille de la matrice soit 41 alors que si q0 < 0, il faut que la taille de la matrice soit 43. Vraiment très bizarre mais j'en ai marre de me prendre la tête sur ça.
        if bec_quasi_momentum < 0:
            size = 43
        else:
            size = 41

        mean_field_energy = (self.atomic_density * self.atom.g / self.E).to(
            u.dimensionless
        )
        V0 = (self.V0 / self.E).to(u.dimensionless)

        def generate_matrix(size, q):
            if size % 2 != 1:
                print("**Warning** Matrix size must be odd")
                return 0
            a = np.zeros((size, size))
            nbsim = int((size - 1) / 2)
            for k in range(size):
                j = k - nbsim
                a[k, k] = (q + 2 * j) ** 2 + (V0 / 2)
                if k != size - 1:
                    a[k, k + 1] = -V0 / 4
                    a[k + 1, k] = -V0 / 4
            #    print(nbsim)
            return a

        def energie(q):
            # q = latt.momentum_to_quasimomentum(q)
            u = generate_matrix(size, q)
            vlp, vpr = LA.eig(u)
            return min(vlp)

        def equations(p):
            q1, q2 = p
            return (
                q1 + q2 - 2 * bec_quasi_momentum + np.sign(bec_quasi_momentum) * 2,
                energie(q1)
                + energie(q2)
                - 2 * energie(bec_quasi_momentum)
                + 2 * mean_field_energy,
            )

        q1, q2 = fsolve(equations, (0.5, 0.5))
        quu = self.momentum_to_quasimomentum(np.array([q1, q2]))
        return (bec_quasi_momentum, np.min(quu), np.max(quu))

    def get_atomic_density(self, z_born=1.5, z_size=1000, bec_speed=0 * u.mm / u.s):
        bec_quasi_momentum = (
            self.atom.speed_to_momentum(bec_speed - self.v) / self.k
        ).to(u.dimensionless)
        Cjmatrix = self._generate_Cjmatrix_dimensionless(bec_quasi_momentum)
        vlp, vpr = LA.eig(
            Cjmatrix
        )  # vlp are eigenvalues --> energies and vpr is the eigenvector associated.
        fondam_cj_coefficients = vpr[:, list(vlp).index(min(vlp))]
        # les coefficients Cj sont rangées de -Cjmatrix_size/2 à Cjmatrix_size/2

        a = pi / 1  # we place our self with dimensionless units therefore 1 = self.k
        z = np.linspace(-z_born * a, z_born * a, z_size)

        def densite(z):
            sum = 0
            for k in range(self.Cjmatrix_size):
                sum = sum + fondam_cj_coefficients[k] * np.exp(
                    1j
                    * z
                    * (bec_quasi_momentum + 2 * np.pi * (k - self.Cjmatrix_size - 1))
                    / a
                )

            return np.abs(sum) ** 2

        densityprofile = densite(z)
        return z, densityprofile

    ######################################################
    ### INTERNAL FUNCTIONS FOR CALCULATIONS
    ######################################################
    def _generate_Cjmatrix_dimensionless(self, q):
        V0 = (self.V0 / self.E).to(u.dimensionless)
        self.check_Cjmatrixsize()

        a = np.zeros((self.Cjmatrix_size, self.Cjmatrix_size))
        nbsim = int((self.Cjmatrix_size - 1) / 2)
        for k in range(self.Cjmatrix_size):
            j = k - nbsim
            a[k, k] = (q + 2 * j) ** 2 + (V0 / 2)
            if k != self.Cjmatrix_size - 1:
                a[k, k + 1] = -V0 / 4
                a[k + 1, k] = -V0 / 4
        return a

    def _generate_Cjmatrix(self, q):
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

    def energy_conservation(self, q1, q2, density=0):

        print(1)
        # g_GP = 1.58e-49  # constante de couplage GP
        # n0 = alpha * Erec / (2 * g_GP)
        return 1

    def get_pairs_quasi_momentum_not_working(
        self,
        atomic_density=1.3e13 / ((1 * u.cm) ** 3),
        bec_speed=0 * u.mm / u.s,
        precision=201,
    ):
        """Cett fonction ne marche pas bien et je ne sais pas pourquoi....

        Parameters
        ----------
        atomic_density : int, optional
            _description_, by default 0
        bec_speed : _type_, optional
            _description_, by default 0*u.mm/u.s
        precision : int, optional
            _description_, by default 201

        Returns
        -------
        _type_
            _description_
        """

        ## D'abord, il faut trouver la relation de disperion de la bande fondamentale
        n_bands = 41  # cf thèse de Maxime
        quasimomentum, energies = self.get_bands_structure(
            n=n_bands, precision=precision
        )
        # /!\/!\/!\
        # Dans la suite de la fonction, on travaille sans dimension !!
        # /!\/!\/!\
        # On traite maintenant le problème SANS DIMENSION car concaténer des arrays avec des dimensions pose problème.
        bec_quasi_momentum = (
            self.atom.speed_to_momentum(bec_speed - self.v) / self.k
        ).to(u.dimensionless)
        # On récupère seulement la bande fondamentale.
        fondam = (energies[0, :] / self.E).to(u.dimensionless)
        fondam = np.concatenate((fondam, fondam, fondam))
        momentum = (quasimomentum / self.k).to(u.dimensionless)
        momentum = np.concatenate((momentum - 2, momentum, momentum + 2))
        ### RECHERCHE DE LA PAIRE
        ## On définit une fonction E (l'énergie du fondamental) qui interpole notre résultat
        E = interpolate.interp1d(momentum, fondam)
        mean_field_energy = (atomic_density * self.atom.g / self.E).to(u.dimensionless)
        # On recherche maintenant les paires (recherche dummy de racine)
        def energy_conservation(q):
            return (
                2 * E(bec_quasi_momentum)
                - E(q)
                - E(2 * bec_quasi_momentum - q)
                + 2 * mean_field_energy
            )

        # Now we want to fin the zero of the function (I cannot use fsolve because the function is not defined evrywhere and it brake the programm.)
        q = np.linspace(bec_quasi_momentum - 1, bec_quasi_momentum + 1, 5 * precision)
        indice_first_zero = np.argmin(
            np.sign(energy_conservation(q[0:-1]) * energy_conservation(q[1:]))
        )
        if indice_first_zero == 0:
            q1 = bec_quasi_momentum
            q2 = bec_quasi_momentum
        else:
            q1 = q[indice_first_zero]
            q2 = 2 * bec_quasi_momentum - q1

        ## CONVERSION MOMENTUM --> QUASIMOMENTUM
        q1 = self.momentum_to_quasimomentum(q1)
        q2 = self.momentum_to_quasimomentum(q2)
        return (bec_quasi_momentum, min(q1, q2), max(q1, q2))


# -- TESTS
if __name__ == "__main__":
    lat = Lattice()
    lat.show_lattice_properties()
