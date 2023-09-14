# !/usr/bin/env python
#  -*-mode:Python; coding:utf8 -*-

# ------------------------------------
# Created on We Sep 2023 by Victor Gondret
# Contact : victor.gondret@institutoptique
#
# MIT Copyright (c) 2023 - Helium1@LCF
# Institut d'Optique Graduate School
# Université Paris-Saclay
# ------------------------------------
#
"""
Decription of parameter_model.py 

Please document your code ;-).
"""
import os, re, logging, json
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


AUTHORIZED_TYPES = [int, float, str, bool, list]


class GlobalModelForParameters:
    def __init__(self, parameters={}, file=__file__):
        self._child_file = file
        self._default_parameters_extension = ".json"
        self._default_parameter_folder = os.path.join(
            os.path.dirname(__file__), "default_parameters"
        )
        self._set_parameter_directory()
        default_parameters = self._load_default_parameters()
        self.__dict__.update(parameters)
        self._save_parameters()

    def get_parameters_str(self):
        """return every attribut of the class as a strig if it is a boolean, a string or a number."""
        parameters = self.get_parameters()
        for key, value in self.__dict__.items():
            key = key.replace("_", " ")
            value = str(value)
        return parameters

    def get_parameters(self):
        parameters = {}
        for key, value in self.__dict__.items():
            if key[0] == "_":
                continue  # les variables commencant par _ sont réservées.
            if type(value) in AUTHORIZED_TYPES:
                parameters[key] = value
        return parameters

    def update_parameter(self, key, str_val):
        """This function update the parameters dictionary, giving a key of the dictionary and a str_val"""
        key = key.replace(" ", "_")
        if key not in self.__dict__.keys():
            logger.warning(f"WARNING : The {key} is not in the model.")
            return str_val
        initial_value = self.__dict__[key]
        try:
            if type(initial_value) == bool:
                if str_val.lower() in "vrai true":
                    self.__dict__[key] = True
                elif str_val.lower() in "faux false":
                    self.__dict__[key] = False
                else:
                    logger.warning(
                        f"WARNING: parameter boolean {key} was not recognized."
                    )
            elif type(initial_value) in [float, int]:
                value = float(str_val)
                if int(value) == value:
                    self.__dict__[key] = int(value)
                else:
                    self.__dict__[key] = value
            elif type(initial_value) == str:
                self.__dict__[key] = str_val
            elif type(initial_value) == list:
                str_val_list = str_val.replace("[", "").replace("]", "").split(",")
                value = []
                for i in range(len(str_val_list)):
                    str_val = str_val_list[i].replace(" ", "")
                    try:
                        val = float(str_val)
                        if int(val) == val:
                            value.append(int(val))
                        else:
                            value.append(val)
                    except:
                        logger.warning(
                            f"Warning: conversion of {str_val} to float in {key} failed."
                        )

                self.__dict__[key] = value
            else:
                logger.warning(
                    f"WARNING : I did not recognized the type of the parameter {key}"
                )
            return str(self.__dict__[key])
        except Exception as e:
            logger.warning(
                f" The entry {key}={str_val} in your parameter is wrong. We ignore this.  \n Error description : {e}"
            )
            return str(self.__dict__[key])

    def _save_parameters(self):
        """
        Sauvegarde les paramètres du modèle dans le fichier de configuration.
        """
        if self._parameter_filename:
            params = self.get_parameters()
            try:
                with open(self._parameter_filename, "w") as file:
                    json.dump(params, file, indent=4)
            except Exception as err:
                logging.error(
                    f" something went wrong trying to save parameters in {self._parameter_filename}"
                )

    def _set_parameter_directory(self):
        current_file_path = os.path.dirname(__file__)
        self._parameter_folder = os.path.join(current_file_path, "parameters")
        if self._child_file == __file__:
            logger.warning(
                "WARNING: the file was not given when you instanciate your model. Please do not forget to give file = __file__ in your super().__init__"
            )
            self._parameter_filename = False
        else:
            # split the child file name as a list
            # par exemple ['', 'home', 'victor', 'heliumtools', 'heliumtools', 'GUIs', 'gaussian_example.py']
            child_file_splitted = os.path.abspath(self._child_file).split(os.sep)
            # we built the name of the default parameter file
            # as GUIs_gaussian_example.json
            self._parameter_filename = ""
            while child_file_splitted:
                prefix = child_file_splitted[-1]
                if prefix in "home heliumtools":
                    break
                self._parameter_filename = prefix + "_" + self._parameter_filename
                child_file_splitted.pop(-1)
            # - supprime l'extension pour mettre la bonne extension
            self._parameter_filename = (
                re.sub(r"\..*", "", self._parameter_filename)
                + self._default_parameters_extension
            )
            self._parameter_filename = os.path.join(
                self._default_parameter_folder, self._parameter_filename
            )

        # file_location = inspect.getfile(my_function)

    def _load_default_parameters(self):
        try:
            if self._parameter_filename:
                with open(self._parameter_filename, "r") as file:
                    default_params = json.load(file)
                logging.info(
                    f"Default parameters loaded from {self._parameter_filename}"
                )
                return default_params
            else:
                return {}
        except Exception as err:
            logging.error(
                f" Error when trying to load parameters from {self._parameter_filename}.  {err}"
            )
            return {}
