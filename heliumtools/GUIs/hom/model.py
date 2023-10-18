# !/usr/bin/env python
#  -*-mode:Python; coding:utf8 -*-

# ------------------------------------
# Created on Mo Sep 2023 by Victor Gondret
# Contact : victor.gondret@institutoptique
#
# MIT Copyright (c) 2023 - Helium1@LCF
# Institut d'Optique Graduate School
# Université Paris-Saclay
# ------------------------------------
#
"""
Decription of model.py 

Model class for the HOM App
"""


import numpy as np

from heliumtools.GUIs.base.parameter_model import GlobalModelForParameters
from heliumtools.misc.hom_database import (
    get_all_BS_delay,
    load_data_BS_delay,
    data_filter,
    load_all_bec_arrival_times,
)
from heliumtools.misc.gather_data import apply_ROI
from math import pi, ceil
import logging
from heliumtools.correlations import Correlation
from tqdm import tqdm

logger = logging.getLogger(__name__)  #


def gaussian(x, mean, amplitude, standard_deviation):
    return amplitude * np.exp(-((x - mean) ** 2) / (2 * standard_deviation**2))


class Model(GlobalModelForParameters):
    """
    The model inherits the GlobalModelForParameters class from heliumtools.misc.
    """

    def __init__(self, hom_folder, parameters={}):
        self._hom_folder = hom_folder
        self.all_delay = get_all_BS_delay(self._hom_folder)
        self.arrival_time = 307.5
        self.inertial_frame = [0, 0, 94]
        self.bragg_speed = 49.6
        self.boxZsize = 1
        self.boxXsize = 15
        self.boxYsize = 15
        self.var_min = 20
        self.var_max = 30
        self.var_step = 0.25
        self.fig_vmax = 1.5
        self.fig_vmin = 0.8
        self.cmap = "viridis"
        self.value_to_plot = "g^2"
        self.Xpos1 = 0
        self.Xpos2 = 0
        self.Ypos1 = 0
        self.Ypos2 = 0
        self.oneD_frame_gap = 5
        self.oneD_curves_per_plot = 4
        self.filter_bec_arrival_time = [307.25, 307.75]
        self.filter_num_of_atoms = [3300, 5500]
        self.filter_num_of_atoms_in_ROI = [100, 280]
        self.filter_m1_atoms = [0, 1000]
        # The parameter dictionary given in the initialisation will overwrite
        # other parameters
        super().__init__(parameters=parameters, file=__file__)
        self._all_bec_arrival_times = (
            load_all_bec_arrival_times(self._hom_folder)
            .sort_values(by="Cycle")
            .reset_index(drop=True)
        )
        self.apply_filters()

    def apply_filters(self):
        """
        Method that update the bec_arrival_times dataframe
        """
        self.filters = {
            "BEC Arrival Time": self.filter_bec_arrival_time,
            "Number of Atoms": self.filter_num_of_atoms,
            "Number of Atoms in ROI": self.filter_num_of_atoms_in_ROI,
            "Number of Atoms in supplementary ROI{}": self.filter_m1_atoms,
        }
        self.bec_arrival_times = apply_ROI(self._all_bec_arrival_times, self.filters)

    def compute_correlations(self):
        self.apply_filters()
        logging.info("Start to compute correlations.")
        ROI = {
            "Vz": {"min": -50, "max": 50},
            "Vy": {"min": -70, "max": 70},
            "Vx": {"max": 25, "min": -75},
        }
        boxes = {
            "1": {
                "Vx": {"size": self.boxXsize, "position": self.Xpos1},
                "Vy": {"size": self.boxYsize, "position": self.Ypos1},
                "Vz": {"size": self.boxZsize, "position": 25},
            },
            "2": {
                "Vx": {"size": self.boxXsize, "position": self.Xpos2},
                "Vy": {"size": self.boxYsize, "position": self.Ypos2},
                "Vz": {"size": self.boxZsize, "position": -25},
            },
        }
        ref_frame = {
            "Vx": self.inertial_frame[0],
            "Vy": self.inertial_frame[1],
            "Vz": self.inertial_frame[2],
        }
        self.result = []
        for i, delay in tqdm(enumerate(self.all_delay), desc="Computing correlations"):
            atoms, bec_arrival_times = load_data_BS_delay(self._hom_folder, delay)
            atoms, bec_arrival_times = data_filter(
                atoms, bec_arrival_times, self.filters
            )

            self.corr = Correlation(
                atoms,
                ROI=ROI,
                boxes=boxes,
                bec_arrival_time=self.arrival_time,
                ref_frame_speed=ref_frame,
            )
            self.corr.define_variable1(
                box="1",
                axe="Vz",
                type="position",
                name="Vz1",
                min=self.var_min,
                max=self.var_max,
                step=self.var_step,
            )
            self.corr.define_variable2(
                box="2",
                axe="Vz",
                type="position",
                name="Vz2",
                values=self.corr.var1.values - self.bragg_speed,
            )
            self.corr.compute_correlations()
            self.corr.result["Beam splitter delay (us)"] = delay
            self.result.append(self.corr.result)
        logging.info("Computation done.")

    def get_properties_for_title(self):
        title = r" $T_{{BEC}}$ = {} ms and $V_{{ref}}$={} mm/s. ".format(
            self.arrival_time, self.inertial_frame
        )
        title += "\n"
        title += r"$\Delta V_z =$ {} mm/s ; $\Delta V_x $= {} mm/s,  $\Delta V_y $= {} mm/s".format(
            self.boxZsize, self.boxXsize, self.boxYsize
        )
        return title
