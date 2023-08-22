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
import textwrap, json
from logging import exception
#from qcontrol3.tools import qc3config


class ConfigFile(object):
    """
    Note : this class is stolen from the qcontrol3.tools.config file.
    A 'normal' usage is to create a config file and to store properties
    in it. We do not do that yet for reasons I have in mind and I am too
    lazy to write.
    
    Base class for a JSON configuration file in the qcontrol3
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

    def __init__(self):
        self.json_data = None
        # The filename stub (without .json extension) for the
        # configuration file represented by this subclass.
        self._filename = "basic_gui"

        # Configuration options and their default values
        self._content = [
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

        # Error messages. Should be overwritten in derived classes. Note
        # that to keep the strings readable, every message here will be
        # passed through textwrap.dedent() and .strip().
        #
        self._messages = {
            "file_not_found_error": """
            I did not find the qcontrol3 configuration file under "{filename}".
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
            raise exception(
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
        print(self._filename)
        #p = get_config_file_path(self._filename)
        #try:
        #     with open(p, "r") as f:
        #         pass
        # except FileNotFoundError:
        #     _message(
        #         self._get_message_string("file_not_found_error").format(
        #             filename=self._filename
        #         )
        #     )
        #     self.write()

        # # Try to parse the file into a temporary variable
        # data = None
        # fail_message = "\n\n" + (
        #     self._get_message_string("fail_message").format(filename=p)
        # )
        # try:
        #     with open(p, "r") as f:
        #         data = json.load(f)
        # except Exception as e:
        #     self._message(self._get_message_string("json_parse_error") + fail_message)
        #     raise e
        data = self.json_data
        #data = self._content
        fail_message = "Something went wrong in base.qt.app.config.py file."
        # Update the self.json_data with data loaded from file, but
        # only if sections and keys are known.

        d_set = set(self.json_data.keys()) - set(data.keys())
        for s_name in d_set:
            self._message(
                self._get_message_string("section_not_found_error").format(
                    section=s_name
                )
                + fail_message
            )

        for s_name, s_items in data.items():
            if s_name not in self.json_data:
                self._message(
                    self._get_message_string("section_unknown_error").format(
                        section=s_name
                    )
                    + fail_message
                )
            else:
                s = self.json_data[s_name]
                d_set = set(s.keys()) - set(s_items.keys())
                for i_name in d_set:
                    self._message(
                        self._get_message_string("option_not_found_error").format(
                            option=i_name, section=s_name, default=s[i_name]
                        )
                        + fail_message
                    )

                for i_name, i_value in s_items.items():
                    if i_name not in s:
                        self._message(
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
        #p = get_config_file_path(self._filename)
        #with open(p, "w") as f:
        #    json.dump(
        #        self.json_data, f, indent=4, separators=(",", " : "), sort_keys=True
        #    )
        pass 

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
            raise exception(
                'Unknown configuration file section "{}". Known sections are {}'.format(
                    section, list(self.json_data)
                )
            )
        if not self.has_option(section, option):
            raise exception(
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
            raise exception(
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
            raise exception(
                'Unknown configuration file section "{}"'.format(section)
            )
        if not self.has_option(section, option):
            raise exception(
                'Unknown configuration file option "{}"'.format(option)
                + ' under section "{}"'.format(section)
            )
        s = self.json_data[section]
        if not isinstance(value, type(s[option])):
            raise exception(
                'Trying to set configuration file option "{}"'.format(option)
                + ' under section "{}" with type "{}"'.format(section, type(s[option]))
                + ' to "{}" with mismatched type "{}"'.format(value, type(value))
            )
        s[option] = value


# config.py ends here
