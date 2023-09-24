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
Decription of example.py 

Example of a heliumtools application. This application use 3 classes : 
    * The HeliumAppBase which is QtWindows whose central widget contains three
    parts. The right one contains tabs with figure, the left top widgets are 
    parameters of the model and the left bottom contains custom buttons 
    defined in the code just above.
    * The model should be defined in a different file, for example here in the
    model.py file. It inherits from the GlobalModelForParameter class that 
    ables the App to get all attributs of the model to show them as parameter.
    * figures that should be shown on the tab widgets inherit from the 
    ClassFigure  class. 

    
Useful method of the HeliumAppBase :     
    * self.add_user_widget(widget)
        add a widget to the left-bottom widget of the app. Give this method 
        any custom button you want to add to your app.
    * self.update_all_plots()
        Use this function to update all the app plots
    * self.update_plot(tab_number)
        To update only one plot defined by its tab_number
    * self.update_model_with_user_parameters()
        This method updates the model with the parameters given by the user
    * self.set_icon(icons = "path_to_your_icon")
        use this method to define your favourite icon as the app icon
"""


import os, logging, sys
import numpy as np
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

# Configuration des logs quand le programme est appelé localement
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)  #


class HOMApp(HeliumAppBase):
    def __init__(self, hom_folder, parameters={}):
        """
        Initialisation of the GUIs. Before initialisation of the GUI, we must create
        """
        # Define the model
        model = Model(hom_folder=hom_folder, parameters=parameters)
        tab_names = [
            "Stability",
            "2D Correlation maps",
            "Hong-Ou-Mandel ?",
        ]
        super().__init__(model=model, tab_names=tab_names, file=__file__)
        # Custom your app
        self.setWindowTitle("Example of an Application.")
        self.setGeometry(100, 100, 1000, 400)
        self.set_icon("gauss.jpeg")
        self.update_correlations_callback()
        self.update_all_plots()

    def setup_user_buttons(self):
        """
        Customize this method to create buttons related to custom callback functions.
        To do so :
            1. Create your button (or widget),
            2. Link it to a callback function that calls a function from your model
            3. Add this button to the 'button pannel' using the add_user_button(button) method
        """
        buttonFilter = QPushButton("Update Parameters")
        buttonFilter.clicked.connect(self.update_filters_callback)
        self.add_user_widget(buttonFilter)

        buttonCorrelations = QPushButton("Compute correlations")
        buttonCorrelations.clicked.connect(self.update_correlations_callback)
        self.add_user_widget(buttonCorrelations)

    ################################################################
    ## ---- Define yere your callbcaks
    ################################################################

    def update_filters_callback(self):
        """We update only filters and the stability plot."""
        self.update_model_with_user_parameters()
        self.model.apply_filters()
        self.update_plot(tab_number=0)

    def update_correlations_callback(self):
        self.update_model_with_user_parameters()
        self.model.compute_correlations()
        self.update_all_plots()


def main(hom_folder, parameters):
    # launch app
    app = QApplication(sys.argv)
    density_app = HOMApp(hom_folder=hom_folder, parameters=parameters)
    density_app.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    HOM_FOLDER = "/mnt/manip_E/2023/09/HOM/02"
    parameters = {}
    main(HOM_FOLDER, parameters)
