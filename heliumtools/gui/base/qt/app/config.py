#!/usr/bin/env python
# -*- mode: Python; coding: utf-8 -*-
# Time-stamp: "2018-03-01 16:31:45 srlab"

#  file       config.py
#  author     Sebastian Blatt
#  created    2018-03-01 13:56:02
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

from heliumtools.misc import heliumtoolsconfig


class ConfigFile(heliumtoolsconfig.ConfigFile):
    _filename = "basic_gui"
    _content = [
        (
            "app",
            "Describes the application",
            [
                ("name", "qcontrol3.gui.qt.app.app", "The application name"),
                ("version", "0.1.0", "The application version"),
                ("description", "This GUI is a demo", "The application description"),
                ("authors", [], "The application authors"),
                ("emails", [], "The application authors' email addresses"),
                ("organization_name", "srlab", "Organization name"),
                ("organization_domain", "srlab.net", "Organization domain"),
            ],
        ),
        (
            "gui",
            "Configures the gui",
            [
                (
                    "geometry",
                    [20, 30, 1024, 768],
                    "Window geometry as [X, Y, WIDTH, HEIGHT]",
                ),
                ("main_window_title", "Basic GUI", "Window title"),
                ("main_page_splitter_sizes", [768, 0], "Main page splitter heights"),
            ],
        ),
    ]


# config.py ends here
