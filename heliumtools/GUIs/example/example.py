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
        self.set_icon("gauss.jpeg")

        
        
    def setup_user_buttons(self):
        """
        Customize this method to create buttons related to custom callback functions. 
        To do so :
            1. Create your button (or widget),
            2. Link it to a callback function that calls a function from your model
            3. Add this button to the 'button pannel' using the add_user_button(button) method
        """
        button = QPushButton('Update My Model')
        button.clicked.connect(self.update_parameters_callback)
        self.add_user_widget(button)
        
        button1 = QPushButton('Update model but only Sinc')
        button1.clicked.connect(self.update_parameters_model_but_only_sinc_callback)
        self.add_user_widget(button1)

    def update_parameters_callback(self):
       self.update_model_with_user_parameters()
       # self.model.compute_stuffs()
       self.update_all_plots()
    
    def update_parameters_model_but_only_sinc_callback(self):
        self.update_model_with_user_parameters()
        self.update_plot(1)
       


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
