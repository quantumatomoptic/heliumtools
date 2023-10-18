#!/usr/bin/env python
# -*- mode: Python; coding: utf-8 -*-
# Time-stamp: "2018-03-01 13:39:55 sb"

#  file       icons.py
#  author     Sebastian Blatt
#  created    2018-03-01 11:21:59
#
#  Copyright (C) 2011 -- 2016 Christoph Gohle, Christian Gross
#  Copyright (C) 2016 -- 2018 Christoph Gohle, Christian Gross, Sebastian Blatt
#
#  This file is part of qcontrol -- The MPQ-developed, python-based
#  experiment control program.
#
#  License:
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#  Commentary:
#

import os

from heliumtools.gui.base.qt import qt


def get_icon_path(subdirectory, filename):
    return os.path.join(os.path.dirname(__file__), subdirectory, filename)


class Icon(qt.QtGui.QIcon):
    def __init__(self, subdirectory, filename):
        path = get_icon_path(subdirectory, filename)
        super().__init__(path)


class Pixmap(qt.QtGui.QPixmap):
    def __init__(self, subdirectory, filename):
        path = get_icon_path(subdirectory, filename)
        super().__init__(path)


class PixmapLabel(qt.QtWidgets.QLabel):
    def __init__(self, parent, subdirectory, filename):
        super().__init__(parent)
        pixmap = Pixmap(subdirectory, filename)
        self.setPixmap(pixmap)


# icons.py ends here
