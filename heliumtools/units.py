#!/usr/bin/env python3
# -*- mode: Python; coding: utf-8 -*-
# Time-stamp: "2018-01-08 18:28:38 srlab"
#
#  file       units.py
#  author     Christoph Gohle
#  created    2012-02-28 03:05:40
#
#  Copyright (C) 2011 -- 2016 Christoph Gohle, Christian Gross
#  Copyright (C) 2016 -- 2019 Christoph Gohle, Christian Gross, Sebastian Blatt
#
#  This file is part of qcontrol3 -- The MPQ-developed, python3-based
#  experiment control program.
#
#  License:
#
#    This program is free software: you can redistribute it and/or
#    modify it under the terms of the GNU General Public License as
#    published by the Free Software Foundation, either version 3 of the
#    License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful, but
#    WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#    General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <http://www.gnu.org/licenses/>.
#

"""
This file provides the units and constanst system for qcontrol3.

"""

from heliumtools import qunits

u = qunits.UnitNamespace("all")
const = qunits.PhysConst(u)

# just for printing we set representations
qunits.setrepresent(u.T, u.G, "G")


from heliumtools.qunits import ScalarQuantity, ArrayQuantity, EIncompatibleUnits

# units.py ends here
