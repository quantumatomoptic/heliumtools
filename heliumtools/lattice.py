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
import numpy as np
from math import pi


class Lattice:
    def __init__(self, **kwargs):
        self.up_frequency = 200 * u.MHz
        self.down_frequency = 200.100 * u.MHz
        self.up_power = 85 * u.mW
        self.down_power = 85 * u.mW
        self.up_waist = 200 * u.um
        self.down_waist = 200 * u.um
        self.theta = 166 / 360 * 2 * pi
        self.wavelength = 1064 * u.nm
        self.__dict__.update(kwargs)

    def build_lattice_properties(self):
        """Construit les propriétés de la classe Lattice à partir des arguments initialisés. Cette fonction est appelée en fin d'initialisation."""
        self.up_intensity = 2 * self.up_power / self.up_waist**2 / pi
        self.down_intensity = 2 * self.down_power / self.down_waist**2 / pi
        self.detuning = self.up_frequency - self.down_frequency
        self.k_latt = np.sin(self.theta / 2) * 2 * pi / self.wavelength
        self.lattice_speed = self.detuning / 2 / self.k_latt
    
    
