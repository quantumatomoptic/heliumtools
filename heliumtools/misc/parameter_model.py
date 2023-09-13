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
import logging
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


class GlobalModelForParameters:
    def __init__(self, parameters={}):
        self.__dict__.update(parameters)

    def get_parameters(self):
        """return every attribut of the class as a strig if it is a boolean, a string or a number."""
        parameters = {}
        for key, value in self.__dict__.items():
            if type(value) in [int, float, str, bool, list]:
                parameters[key.replace("_", " ")] = str(value)
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
                str_val = str_val.replace("[", "").replace("]", "").split(",")
                value = []
                for i in str_val:
                    try:
                        if int(i) == float(i):
                            value.append(int(i))
                        else:
                            value.append(float(i))
                    except:
                        logger.warning(
                            f"Warning: conversion of {i} to float in {key} failed."
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
