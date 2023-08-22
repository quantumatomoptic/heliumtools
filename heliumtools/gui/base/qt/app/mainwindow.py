#!/usr/bin/env python
# -*- mode: Python; coding: utf-8 -*-
# Time-stamp: "2018-07-31 11:22:49 srlab"

#  file       mainwindow.py
#  author     Sebastian Blatt
#  created    2018-03-01 13:59:28
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

import sys

from heliumtools.gui.base.qt import qt
from heliumtools.gui.base.qt.widget.dialog import about


# Handling Ctrl-C (SIGINT) cleanly is not easy. See e.g. the
# discussion on
#
#   https://stackoverflow.com/questions/4938723/what-is-the-correct-way-to-make-my-pyqt-application-quit-when-killed-from-the-co
#   https://coldfix.de/2016/11/08/pyqt-boilerplate/
#
# The problem with the code shown there is that a QTimer can only be
# started within a QThread. This means that we have to start the
# safe_timer function from within the main window.
#
# In addition, we want to be able to gracefully close the window
# instead of simply killing the app. This means that sigint_handler
# needs access to the MainWindow instance.
#


def safe_timer(timeout, func, *args, **kwargs):
    def timer_event():
        try:
            func(*args, **kwargs)
        finally:
            qt.QtCore.QTimer.singleShot(timeout, timer_event)

    qt.QtCore.QTimer.singleShot(timeout, timer_event)


class MainWindow(qt.QtWidgets.QMainWindow):
    aboutToCloseEvent = qt.QtCore.pyqtSignal()

    def __init__(self, model, main_widget_type):
        super().__init__()
        self._model = model
        self._main_widget_type = main_widget_type
        self._main_widget = None
        self._actions = {}
        self._initialize_ui()

        # This magic timer is used to handel SIGINT cleanly. See
        # below.
        safe_timer(200, lambda: None)

        self._model.start_threads()

    def _initialize_ui(self):
        self._setup_main_widget()
        self._setup_actions()

        self.statusBar()
        self._setup_menubar()

        self.setGeometry(*self._model.get_config("gui", "geometry"))
        self.setWindowTitle(self._model.get_config("gui", "main_window_title"))
        self.show()

    def _setup_main_widget(self):
        self._main_widget = self._main_widget_type(self, self._model)
        self.setCentralWidget(self._main_widget)

    def _setup_actions(self):
        action_exit = qt.QtWidgets.QAction("Quit", self)
        action_exit.setShortcut("Ctrl+Q")
        action_exit.setStatusTip("Exit Application")
        action_exit.triggered.connect(self.close)
        self._actions["exit"] = action_exit

        action_about = qt.QtWidgets.QAction("About", self)
        action_about.setStatusTip(
            "About {}".format(self._model.get_config("app", "name"))
        )
        action_about.triggered.connect(self._on_about)
        self._actions["about"] = action_about

    def _setup_menubar(self):
        if sys.platform == "darwin":
            self.menuBar().setNativeMenuBar(False)
        self._setup_menu_file()
        self._setup_menu_help()

    def _setup_menu_file(self):
        m = self.menuBar()
        fm = m.addMenu("&File")
        fm.addAction(self._actions["exit"])
        return fm

    def _setup_menu_help(self):
        m = self.menuBar()
        hm = m.addMenu("&Help")
        hm.addAction(self._actions["about"])
        return hm

    def _on_about(self):
        dialog = about.AboutDialog(self, self._model)
        dialog.setWindowModality(qt.QtCore.Qt.ApplicationModal)
        dialog.exec_()

    def closeEvent(self, *args, **kwargs):
        super().closeEvent(*args, **kwargs)
        x = self.geometry()
        self._model.set_config(
            "gui", "geometry", [x.left(), x.top(), x.width(), x.height()]
        )

        # hide main window immediately to let the user know that
        # "something is being done"
        self.hide()

        # handle GUI shutdown gracefully
        self.aboutToCloseEvent.emit()

        self._model.stop_threads()
        self._model.write_config()


# mainwindow.py ends here
