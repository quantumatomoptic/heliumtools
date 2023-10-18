#!/usr/bin/env python
# -*- mode: Python; coding: utf-8 -*-
# Time-stamp: "2018-08-03 13:06:27 srlab"

#  file       model.py
#  author     Sebastian Blatt
#  created    2018-03-01 13:57:03
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
from heliumtools.gui.base.qt.widget.logging import logtable
from datetime import datetime


class Model:
    def __init__(self, config_file, log_storage=None):
        # we will store the global application object and the main window here
        self._application = None
        self._main_window = None

        # the model manages the configuration file interface
        self._config_file = config_file
        self._config_file.read()

        # as well as the internal log storage
        if log_storage is None:
            log_storage = logtable.LogStorage(256)

        self._log_storage = log_storage

        # the model manages all daughter threads
        self._threads = []
        self._thread_stop_timeout_ms = 2000

    @property
    def application(self):
        return self._application

    @application.setter
    def application(self, value):
        if self._application is not None:
            raise Exception
        self._application = value

    @property
    def main_window(self):
        return self._main_window

    @main_window.setter
    def main_window(self, value):
        if self._main_window is not None:
            raise Exception
        self._main_window = value

    def get_config(self, section, option):
        return self._config_file.get(section, option)

    def set_config(self, section, option, value):
        self._config_file.set(section, option, value)

    def write_config(self):
        self._config_file.write()

    def log(self, message, level="DEBUG"):
        self._log_storage.append(
            {
                "levelname": level,
                "asctime": datetime.now().strftime("%H:%M:%S"),
                # 'pathname': '',
                # 'funcName': '',
                # 'lineno': 0,
                # 'process': '',
                # 'threadName': '',
                "msg": message,
            }
        )

    def start_threads(self):
        for t in self._threads:
            t.start()

    def stop_threads(self):
        for t in self._threads:
            t.stop()
        for t in self._threads:
            t.wait(self._thread_stop_timeout_ms)
        for t in self._threads:
            t.quit()


# model.py ends here
