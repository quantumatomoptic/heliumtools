#!/usr/bin/env python
# -*- mode: Python; coding: utf-8 -*-
# Time-stamp: "2018-02-26 17:12:47 srlab"

#  file       plotcanvas.py
#  author     Sebastian Blatt
#  created    2018-02-26 15:41:05
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

from heliumtools.gui.base.qt import qt

from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)

from matplotlib.figure import Figure
from matplotlib import pyplot as plt


class PlotWidget(FigureCanvas):
    aboutToCloseEvent = qt.QtCore.pyqtSignal()

    def __init__(self, parent, model, width=4, height=4, dpi=100):

        self._fig = Figure(figsize=(width, height), dpi=dpi)

        super().__init__(self._fig)
        self.setParent(parent)

        self.mpl_connect('button_press_event', self._onclick)

        self._model = model
        self._initialize_ui()
        parent.aboutToCloseEvent.connect(self._on_about_to_close)

    def _initialize_ui(self):
        self.setSizePolicy(
            qt.QtWidgets.QSizePolicy.Expanding, qt.QtWidgets.QSizePolicy.Expanding
        )
        self.updateGeometry()
        self.plot()
        self.show()

    def _onclick(self, event):
        if event.xdata is not None and event.ydata is not None:
            print('x=%d, y=%d' % (event.xdata, event.ydata))

    def _on_about_to_close(self):
        self._save_geometry()
        self.aboutToCloseEvent.emit()

    def plot(self, data=None):
        pass

    def _save_geometry(self):
        pass


# plotcanvas.py ends here
