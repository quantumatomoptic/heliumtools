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

Please document your code ;-).

"""
import scipy.constants as cst

# from qcontrol3.tools.units import u


class Lattice:
    def __init__(self, **kwargs):
        self.up_frequency = 1
        self.down_frequency = 1
        self.up_power = 1
        self.__dict__.update(kwargs)
