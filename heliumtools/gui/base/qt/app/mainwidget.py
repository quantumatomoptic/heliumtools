#!/usr/bin/env python
# -*- mode: Python; coding: utf-8 -*-
# Time-stamp: "2018-03-01 16:25:49 srlab"

#  file       mainwidget.py
#  author     Sebastian Blatt
#  created    2018-03-01 13:58:31
#
#  Copyright (C) 2011 -- 2016 Christoph Gohle, Christian Gross
#  Copyright (C) 2016 -- 2018 Christoph Gohle, Christian Gross, Sebastian Blatt
#
#  This file is part of qcontrol3 -- The MPQ-developed, python3-based
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

from . import page


class MainWidget(qt.QtWidgets.QTabWidget):
    aboutToCloseEvent = qt.QtCore.pyqtSignal()

    def __init__(self, parent, model):
        super().__init__(parent)
        self._model = model
        self._pages = []
        self._initialize_ui()
        parent.aboutToCloseEvent.connect(self._on_about_to_close)

    def _initialize_ui(self):
        self._pages = self._create_pages()
        for page in self._pages:
            self.addTab(page, page.name)

    def _create_pages(self):
        return [
            page.MainPage(self, self._model),
            page.PreferencesPage(self, self._model),
        ]

    def _on_about_to_close(self):
        self._save_geometry()
        self.aboutToCloseEvent.emit()

    def _save_geometry(self):
        pass


# mainwidget.py ends here
