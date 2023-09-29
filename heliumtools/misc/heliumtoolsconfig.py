#!/usr/bin/env python3
# -*- mode: Python; coding: utf-8 -*-
# Time-stamp: "2018-10-08 15:49:13 srlab"
#
#  file       heliumtoolsconfig.py
#  author     Sebastian Blatt
#  created    2017-02-21 13:44:09
#
#  Copyright (C) 2011 -- 2016 Christoph Gohle, Christian Gross
#  Copyright (C) 2016 -- 2019 Christoph Gohle, Christian Gross, Sebastian Blatt
#
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
Configuration infrastructure for heliumtools.

We would like to store all heliumtools-related configuration files
in the same location, some version of the posix directory ~/.heliumtools,
but portable across platforms. We would like to have independent
configuration files for different heliumtools components.

This file implements infrastructure to create, read, and write
such configuration file in JSON format.

The default directory for storing the JSON files is

  ~/.heliumtools

where the tilde is expanded using $HOME on the windows platform
as well and the leading dot is removed. If $HOME is undefined
under windows, the tilde is expanded to "C:". Under windows, the
other heliumtools directories are *NOT* put in separate directories
under HOME, but are placed underneath the qc3 directory for
better compatibility with windows conventions.

"""

import os
import sys
import shutil
import textwrap
import json

from . import exceptions


# String identifiers for supported platforms
CONFIG_SUPPORTED_PLATFORMS = ("linux", "osx", "windows")

CONFIG_DIRECTORY_BASENAME = ".heliumtools_configuration"

CONFIG_QC3_CONFIG_FILE = "heliumtoolsconf"


# Try to do highly visible logging via sys.stderr and sys.stdout since
# we cannot rely on logging facility yet (it is being configured here).


def _message_header(string=None):
    """Return a line of dashes as wide as the current terminal optionally
    containing the header `string` after the first three dashes.
    """
    dim = shutil.get_terminal_size()
    w = dim.columns

    if string is None or len(string) == 0:
        return "-" * w

    t = "--- " + string
    t += " " + ("-" * max(0, w - len(t) - 1))
    return t + "\n"


def _message_footer():
    """Return a line of dashed as wide as the current terminal."""
    dim = shutil.get_terminal_size()
    return ("-" * dim.columns) + "\n"


def _message(msg, header="Heliumtools initial configuration message"):
    """Wrap message `msg` using _message_header() and _message_footer()
    and print to sys.stdout."""
    print(_message_header(header) + msg + "\n" + _message_footer(), file=sys.stdout)


def _error(msg, header="qcontrol3 initial configuration error"):
    """Wrap error `msg` using _message_header() and _message_footer() and
    print to sys.stderr."""
    print(_message_header(header) + msg + "\n" + _message_footer(), file=sys.stderr)


def get_platform():
    """Determine the OS we are running on by looking at os.name and
    sys.platform. Return string 'linux', 'darwin', or 'windows',
    otherwise raise Exception to indicate an unsupported platform.
    """
    if os.name == "posix":
        if sys.platform == "linux" or sys.platform == "linux2":
            return "linux"
        elif sys.platform == "darwin":
            return "osx"
    elif os.name == "nt" and sys.platform == "win32":
        return "windows"
    else:
        raise exceptions.ConfigFileException(
            "Unknown operating system (os.name, sys.platform) = ({}, {}).".format(
                os.name, sys.platform
            )
        )


def get_absolute_path(*args):
    """Handle all path name conversions in the same way: assume that we
    get a path as a string or as separate arguments. If separate, join
    them using os.path.join. Then expand user directories using
    os.path.expanduser and normalize the path using os.path.abspath.
    """
    p = None
    if len(args) > 1:
        p = os.path.join(*args)
    else:
        p = args[0]
    return os.path.abspath(os.path.expanduser(p))


CONFIG_DIRECTORY_PLATFORM_NAMES = {
    "linux": get_absolute_path("~", "." + CONFIG_DIRECTORY_BASENAME),
    "osx": get_absolute_path("~", "." + CONFIG_DIRECTORY_BASENAME),
    "windows": get_absolute_path("C:", CONFIG_DIRECTORY_BASENAME),
}


def get_config_directory_path():
    """The configuration directory is defined in a platform-independent
    way. For a linux or osx system, we will put the directory in the
    home directory with a leading dot to indicate a hidden directory.
    However, for windows systems, we do not want to assume that $HOME
    is defined and will put the directory under "C:" without a leading
    dot. But, if $HOME is defined, we place the all qc3 directories
    there instead of under "C:".
    """
    p = get_platform()
    path = CONFIG_DIRECTORY_PLATFORM_NAMES[p]
    print(path)
    if (p == "windows") and (os.getenv("HOME") is not None):
        path = get_absolute_path(os.getenv("HOME"), CONFIG_DIRECTORY_BASENAME)
    return path


def get_config_file_path(filename):
    """Return an absolute path for `filename` assuming it represents a
    filename (without file extension) in the heliumtools configuration
    directory.
    Example:
      get_config_file_path('myclient') -> '/home/myuser/.qc3/myclient.json'
    Raises a ConfigFileException if the filename is specified with an extension.
    """
    root, ext = os.path.splitext(filename)
    if ext != "":
        raise exceptions.ConfigFileException(
            "Configuration filenames must be specified without "
            'file extension, i.e. "myfile" instead of "myfile.json"'
        )
    

    cdir = get_config_directory_path()
    p = get_absolute_path(cdir, filename + ".json")
    return p


def check_config_file_permissions(filename):
    """Check that `filename` is a readable and writeable regular file
    within the heliumtools configuration directory or raise a
    ConfigFileException.
    """
    p = get_config_file_path(filename)
    if not os.path.exists(p):
        raise exceptions.ConfigFileException('Path "{}" does not exist'.format(p))
    if not os.path.isfile(p):
        raise exceptions.ConfigFileException(
            'Path "{}" is not a regular file'.format(p)
        )
    if not os.access(p, os.R_OK):
        raise exceptions.ConfigFileException('Path "{}" is not readable'.format(p))
    if not os.access(p, os.W_OK):
        raise exceptions.ConfigFileException('Path "{}" is not writeable'.format(p))


def make_directory_path(*args, **kwargs):
    """Create a directory and, if they do not exist already, its parent
    directories, specified as positional arguments (*args). Then check
    that the directory was successfully created and check that we have
    the correct access permissions using the keyword arguments
    (**kwargs).

    Keyword arguments:

    - read = True: If True, check that os.access using os.R_OK
      returns True.

    - write = True: If True, check that os.access using os.W_OK
      returns True

    Example:

       make_directory_path("~", "test", "directory",
                           read=True, write=True)

    Return:

       None, but raises Exception if any requested operation failed.

    """
    p = get_absolute_path(*args)
    if not os.path.exists(p):
        _message('Directory "{}" does not exist yet, I will create it."'.format(p))
        os.makedirs(p)
        if not os.path.exists(p):
            raise exceptions.ConfigFileException(
                'Failed creating directory "{}"'.format(p)
            )

    if not os.path.isdir(p):
        raise exceptions.ConfigFileException(
            'Path "{}" exists but is not a directory'.format(p)
        )

    if "read" in kwargs:
        a = kwargs.pop("read")
        if a and not os.access(p, os.R_OK):
            raise exceptions.ConfigFileException(
                'Do not have read access to directory "{}"'.format(p)
            )

    if "write" in kwargs:
        a = kwargs.pop("write")
        if a and not os.access(p, os.W_OK):
            raise exceptions.ConfigFileException(
                'Do not have write access to directory "{}"'.format(p)
            )


class ConfigFile(object):
    """Base class for a JSON configuration file in the heliumtools
    configuration directory with default values and comments. Since
    JSON does not support a comment syntax, we auto-generate special
    string variables by appending '_comment' to each regular
    configuration option that also contains the default setting for
    the commented variable.

    The ConfigFile object methods are modelled after the
    configparser.ConfigParser interface, but we only implement the
    functionality that we require.

    ConfigFile *MUST* be subclassed and the class variables _filename
    and _content must be specified.

    """

    # The filename stub (without .json extension) for the
    # configuration file represented by this subclass.
    _filename = None

    # Configuration options and their default values
    _content = []

    # Error messages. Should be overwritten in derived classes. Note
    # that to keep the strings readable, every message here will be
    # passed through textwrap.dedent() and .strip().
    #
    _messages = {
        "file_not_found_error": """
        I did not find the heliumtools configuration file under "{filename}".
        'I will create it using default values, but you must modify these
        'according to your setup. See the documentation under "doc/configuration"
        'on how to do this.""",
        "fail_message": """
        There was something wrong with your configuration file "{filename}".
        Please read up on the file's format in the documentation under
        "doc/configuration", fix your configuration file, and try again.
        """,
        "json_parse_error": """
        I failed to parse the file.
        """,
        "section_not_found_error": """
        The configuration section "{section}" is not set in your file.
        I will create it with default values and move on.
        """,
        "section_unknown_error": """
        You have an unknown section "{section}" in your file, I will ignore it and move on.
        """,
        "option_not_found_error": """
        The configuration option "{option}" under section "{section}" is not set in your file.
        I will create it with its default value "{default}" and move on.
        """,
        "option_unknown_error": """
        You have an unknown option "{option}" under section "{section}" in your file.
        I will ignore it and move on.
        """,
    }

    def __init__(self):
        self.json_data = None
        self._populate_json_defaults()
        self.read()

    def _get_message_string(self, identifier):
        return textwrap.dedent(self._messages[identifier]).strip()

    def _populate_json_defaults(self):
        """Parse the class member `_content` and generate a commented
        JSON structure with default values from it.
        """
        self.json_data = {}
        if self._content is None or len(self._content) == 0:
            raise exceptions.ConfigFileException(
                "Programmer brain damage in ConfigFile subclass. "
                + "Must have _content class member that is not None."
            )

        for (s_name, s_comment, s_items) in self._content:
            section = {"_comment": s_comment}
            for i_name, i_default, i_comment in s_items:
                section[i_name] = i_default
                json_default = json.dumps(i_default)
                section[
                    i_name + "_comment"
                ] = i_comment + " [Default value is {}]".format(json_default)
            self.json_data[s_name] = section

    def read(self):
        """Parse the JSON configuration file specified by self._filename and
        provide explicit error messages if something goes wrong.
        """
        p = get_config_file_path(self._filename)
        try:
            with open(p, "r") as f:
                pass
        except FileNotFoundError:
            _message(
                self._get_message_string("file_not_found_error").format(
                    filename=self._filename
                )
            )
            self.write()

        # Try to parse the file into a temporary variable
        data = None
        fail_message = "\n\n" + (
            self._get_message_string("fail_message").format(filename=p)
        )
        try:
            with open(p, "r") as f:
                data = json.load(f)
        except Exception as e:
            _message(self._get_message_string("json_parse_error") + fail_message)
            raise e

        # Update the self.json_data with data loaded from file, but
        # only if sections and keys are known.
        d_set = set(self.json_data.keys()) - set(data.keys())
        for s_name in d_set:
            _message(
                self._get_message_string("section_not_found_error").format(
                    section=s_name
                )
                + fail_message
            )

        for s_name, s_items in data.items():
            if s_name not in self.json_data:
                _message(
                    self._get_message_string("section_unknown_error").format(
                        section=s_name
                    )
                    + fail_message
                )
            else:
                s = self.json_data[s_name]
                d_set = set(s.keys()) - set(s_items.keys())
                for i_name in d_set:
                    _message(
                        self._get_message_string("option_not_found_error").format(
                            option=i_name, section=s_name, default=s[i_name]
                        )
                        + fail_message
                    )

                for i_name, i_value in s_items.items():
                    if i_name not in s:
                        _message(
                            self._get_message_string("option_unknown_error").format(
                                option=i_name, section=s_name
                            )
                            + fail_message
                        )
                    elif i_name[-8:] != "_comment":
                        self.set(s_name, i_name, i_value)

    def write(self):
        """Write the current configuration in JSON format to the configuration
        file specified by self._filename.
        """
        p = get_config_file_path(self._filename)
        with open(p, "w") as f:
            json.dump(
                self.json_data, f, indent=4, separators=(",", " : "), sort_keys=True
            )

    def __repr__(self):
        return json.dumps(
            self.json_data, indent=4, separators=(",", " : "), sort_keys=True
        )

    def has_section(self, section):
        """Check whether `section` is a valid configuration file section."""
        if self.json_data is not None:
            return section in self.json_data

    def has_option(self, section, option):
        """Check whether (`section`, `option`) represents a valid
        configuration file section and option."""
        if self.has_section(section):
            s = self.json_data[section]
            return option in s
        return False

    def get(self, section, option):
        """Return the value of the configuration file option specified by
        (`section`, `option`). Raise a ConfigFileException if either
        `section` or `option` is unknown.
        """
        if not self.has_section(section):
            raise exceptions.ConfigFileException(
                'Unknown configuration file section "{}". Known sections are {}'.format(
                    section, list(self.json_data)
                )
            )
        if not self.has_option(section, option):
            raise exceptions.ConfigFileException(
                'Unknown configuration file option "{}"'.format(option)
                + ' under section "{}"'.format(section)
            )
        s = self.json_data[section]
        return s[option]

    def getboolean(self, section, option):
        """Return the value of the configuration file option specified by
        (`section`, `option`) as boolean. Raise a ConfigFileException if either
        `section` or `option` is unknown or if value cannot be converted to boolean.
        """
        booltrue = set((1, "1", "True", "true"))
        boolfalse = set((0, "0", "False", "false"))
        strval = self.get(section, option)
        if strval in booltrue:
            return True
        elif strval in boolfalse:
            return False
        else:
            raise exceptions.ConfigFileException(
                'Boolean equivalent of {} (option "{}"'.format(strval, option)
                + ' in section "{}")'.format(section)
            )

    def set(self, section, option, value):
        """Set the value of the configuration file option specified by
        (`section`, `option`). Raise a ConfigFileException if either
        `section` or `option` is unknown or if the type of option does
        not match the type we already have (set via defaults).

        """
        if not self.has_section(section):
            raise exceptions.ConfigFileException(
                'Unknown configuration file section "{}"'.format(section)
            )
        if not self.has_option(section, option):
            raise exceptions.ConfigFileException(
                'Unknown configuration file option "{}"'.format(option)
                + ' under section "{}"'.format(section)
            )
        s = self.json_data[section]
        if not isinstance(value, type(s[option])):
            raise exceptions.ConfigFileException(
                'Trying to set configuration file option "{}"'.format(option)
                + ' under section "{}" with type "{}"'.format(section, type(s[option]))
                + ' to "{}" with mismatched type "{}"'.format(value, type(value))
            )
        s[option] = value

def ensure_config_directory_exists():
    """Ensure that the configuration directory returned by
    get_config_directory_path() exists and that we have read and write
    access.

    Because the config module also provides access functions to the
    configuration directory, this function is called whenever the
    module is loaded. See bottom of config.py for how this is done.
    """
    print()
    make_directory_path(get_config_directory_path(), read=True, write=True)


# Stuff that happens on import of this module.
ensure_config_directory_exists()

if __name__ == "__main__":
    pass

# heliumtoolsconfig.py ends here
