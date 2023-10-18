#!/usr/bin/env python3
# -*- mode: Python; coding: utf-8 -*-
# Time-stamp: "2018-08-01 12:29:37 srlab"

#  file       app.py
#  author     Sebastian Blatt
#  created    2018-03-01 10:59:47
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

import sys
import signal

import qdarkstyle

from heliumtools.gui.base.qt import qt
from heliumtools.gui.base.qt.icons import icons
from heliumtools.gui.base.qt.widget.logging import logtable

from . import config
from . import model
from . import mainwidget
from . import mainwindow


class Qcontrol3App:
    def __init__(self, model, main_window_type, main_widget_type):
        self._setup_signal_handlers()

        a = qt.QtWidgets.QApplication(sys.argv)
        icon = icons.Icon("qcontrol3_logo", "256.png")
        a.setWindowIcon(icon)
        a.setStyleSheet(qdarkstyle.load_stylesheet_pyqt5())
        a.setOrganizationName(model.get_config("app", "organization_name"))
        a.setOrganizationDomain(model.get_config("app", "organization_domain"))
        a.setApplicationName(model.get_config("app", "name"))

        self._model = model
        self._model.application = a
        self._main_window = main_window_type(model, main_widget_type)
        self._model.main_window = self._main_window

    def _setup_signal_handlers(self):
        def handle_sigint(signal, frame):
            self._on_sigint(signal, frame)

        signal.signal(signal.SIGINT, handle_sigint)

    def _on_sigint(self, signal, frame):
        sys.stderr.write("\r")
        self._main_window.close()

    def run(self):
        sys.exit(self._model.application.exec_())


if __name__ == "__main__":
    config_file = config.ConfigFile()
    log_storage = logtable.LogStorage(256)
    model = model.Model(config_file, log_storage)
    app = Qcontrol3App(model, mainwindow.MainWindow, mainwidget.MainWidget)
    app.run()


# app.py ends here
