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
Decription of gaussian_example.py 

"""


import os, logging, sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from PyQt5.QtWidgets import (
    QApplication,
    QWidget,
    QLabel,
    QLineEdit,
    QPushButton,
    QVBoxLayout,
    QHBoxLayout,
    QSizePolicy,
    QScrollArea,
    QSplitter,
)
from heliumtools.GUIs.base.app_base import HeliumAppBase
from heliumtools.misc.gather_data import apply_ROI
from model import Model

# Configuration des logs quand le programme est appelé localment
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class ExampleApplication(HeliumAppBase):
    def __init__(self, data=[], metadata={}, parameters={}):
        """
        Initialisation of the GUIs. Before initialisation of the GUI, we must create
        """
        # Define the model
        model = Model(data=data, metadata=metadata, parameters=parameters)
        tab_names = ["Gaussian Curve", "Sinc Curve", "Thomas Fermi profile", "Failed"]
        super().__init__(model=model, tab_names=tab_names, file=__file__)

        self.setWindowTitle("Example of an Application.")
        self.setGeometry(100, 100, 1000, 400)
        self.set_icon("mon_image.png")

        self.initUI()
        # connect the FiguresClass

    def initUI(self):
        """
        Main function that creates the User Interface.
        """
        self._setup_UI_left_part()
        self._setup_UI_right_part()
        self._main_widget = QSplitter()
        self._main_widget.addWidget(self.left_widget)
        self._main_widget.addWidget(self.right_widget)
        self.setCentralWidget(self._main_widget)

    def _setup_UI_left_part(self):
        """This method sets up the left part of the UI. Do not forget to include the parameter_widget_wrapper"""
        self.left_widget = self.parameter_widget

    def _setup_UI_right_part(self):
        """
        sets up the right pannel of your GUI
        """
        # Right part of the windows
        self.right_widget = self.tabs_fig_widget


def main(data, metadata, parameters):
    # launch app
    app = QApplication(sys.argv)
    density_app = ExampleApplication(
        data=data, metadata=metadata, parameters=parameters
    )
    density_app.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    data = []
    metadata = []
    parameters = {}
    main(data, metadata, parameters)
