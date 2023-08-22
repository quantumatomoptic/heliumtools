#!/usr/bin/env python
# -*- mode: Python; coding: utf-8 -*-
# Time-stamp: "2018-10-08 15:48:12 srlab"
#
#  file       exceptions.py
#  author     Christian Gross
#  created    2017-10-02
#
#  Copyright (C) 2011 -- 2017 Christoph Gohle, Christian Gross
#  Copyright (C) 2016, 2017 Christoph Gohle, Christian Gross, Sebastian Blatt
#
#  This file is part of qcontrol -- The MPQ-developed, python3-based
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

"""The exception module contains custom exceptions used in qcontrol3.

"""

import os.path as osp
import sys


class QcontrolException(Exception):
    def __init__(self, *args, **kwargs):
        # make a beep
        sys.stdout.write("\a")
        super().__init__(*args, **kwargs)

    def formatted_message(self):
        raise NotImplementedError


class ConfigFileException(QcontrolException):
    """Exception class thrown by the qcontrol3 configuration system.

    """

    def formatted_message(self):
        pass


class ConfigurationError(QcontrolException):
    """Exception thrown for system configuration errors.

    Use the description to explain what happend.
    """

    def __init__(self, description):
        super().__init__(description)
        self._description = description

    def formatted_message(self):
        """A formatted error message that is intended to be presented to the user.

        """
        msg = "ConfigurationError: " + self._description
        return msg


class VerificationError(QcontrolException):
    """Exception thrown when channel constraints are violated.

    Use the description to explain what happend.
    """

    def __init__(self, description):
        super().__init__(description)
        self._description = description

    def formatted_message(self):
        """A formatted error message that is intended to be presented to the user.

        """
        msg = "VerificationError: " + self._description
        return msg


class ScriptCollisionError(QcontrolException):
    """Exception thrown in case of a channel collision

    Parameters
    ----------
    tracebackinfo : dict
        A dictionary with a list of collision pairs (again a list with three
        elements). The first element is the channel name, the second and third
        the actual collision pair. Each colliding event is described by a
        dictionary.
        tracebackinfo[channel][0][0]['file'] : the source file
        tracebackinfo[channel][0][0]['line'] : the line number
        tracebackinfo[channel][0][0]['command'] : the command in the source

    """

    def __init__(self, tracebackinfo):
        super().__init__()
        self._tracebackinfo = tracebackinfo

    def formatted_message(self):
        """A formatted error message that is intended to be presented to the user.

        """
        msg = "ScriptCollisionError: " + "\n"
        for node_name in self._tracebackinfo.keys():
            for tb in self._tracebackinfo[node_name]:
                msg += "Collision below node {}:\n".format(node_name)
                msg += "Channels {} and {} collide:\n".format(tb[0], tb[2])
                msg += "File {} in line {}: {}\nFile {} in line {}: {}\n".format(
                    osp.basename(tb[1]["file"]),
                    tb[1]["line"],
                    tb[1]["command"].strip(),
                    osp.basename(tb[3]["file"]),
                    tb[3]["line"],
                    tb[3]["command"].strip(),
                )
        return msg


class ScriptError(QcontrolException):
    """Exception thrown for any common script error.

    Use the description to explain what happend.
    """

    def __init__(self, description, tracebackinfo):
        super().__init__(description)
        self._description = description
        self._tracebackinfo = tracebackinfo

    def formatted_message(self):
        """A formatted error message that is intended to be presented to the user.

        """
        msg = "ScriptError: " + self._description + "\n"
        err_no = 0
        for tb in self._tracebackinfo:
            err_no += 1
            msg += "Error number {}:\n".format(err_no)
            msg += "File {} in line {}: {}\n".format(
                osp.basename(tb["file"]), tb["line"], tb["command"].strip()
            )
        return msg


class ScriptObjectError(QcontrolException):
    """Exception to be thrown in a ScriptObject

    Use the description to explain what happend.
    """

    def __init__(self, description):
        super().__init__(description)
        self._description = description

    def formatted_message(self):
        """A formatted error message that is intended to be presented to the user.

        """
        msg = "ScriptObjectError: " + self._description
        return msg


# exceptions.py ends here
