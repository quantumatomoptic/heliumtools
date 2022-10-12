#!/usr/bin/env python
# -*- mode:Python; coding: utf-8 -*-

# ----------------------------------
# Created on the Tue Oct 11 2022 by Victor
#
# Developped by Victor, ...
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
    k_latt : vecteur d'onde du réseau (qunits)
    lattice_speed : vitesse du réseau (qunits)
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
        self.__dict__.update(kwargs)
        self.build_lattice_properties()

    def build_lattice_properties(self):
        """Construit les propriétés de la classe Lattice à partir des arguments initialisés. Cette fonction est appelée en fin d'initialisation."""
        self.up_intensity = 2 * self.up_power / self.up_waist**2 / pi
        self.down_intensity = 2 * self.down_power / self.down_waist**2 / pi
        self.detuning = self.up_frequency - self.down_frequency
        self.k_latt = np.sin(self.theta / 2) * 2 * pi / self.wavelength
        self.lattice_speed = self.detuning / 2 / self.k_latt
        if self.down_intensity != self.up_intensity:
            warnings.warn(
                "Up and down arms do not have the same intensity. I will take the mean to define the intensity. "
            )
        self.intensity = np.sqrt(self.down_intensity * self.up_intensity)
        self.V0 = const.hbar * self.atom.atomic_width**2 * self.intensity

    def get(self):
        return self.__dict__
