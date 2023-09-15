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



The principle is he following : when the application is run, the model Model is instanciated. 
Once it is instanciate, the application gets back all parameters of the model from the function self.model.get_parameters_str(). 
This generates a dictionary from model.__dict__ in which each element is a number/list/string or boolean. 
From this dictionary, the application built in a scrollable area a list of QLabels and Qlines to updtate parameters of the model. 
The figure that is shown on the right is defined in the PlotZaxis class. It does not contains a lot but the method update_plot that is called each time the user pushes the button 'Update Plot'.  This function obviously need the model to be rightly updated.
"""


import os, logging, sys
import numpy as np
import matplotlib.pyplot as plt
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
from PyQt5.QtGui import QIcon, QPixmap, QImage, QClipboard
from PyQt5.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import pandas as pd
from heliumtools.correlations import Correlation
from scipy.optimize import curve_fit
import seaborn as sns
from PIL import Image
from heliumtools.GUIs.base.parameter_model import GlobalModelForParameters
from heliumtools.GUIs.base.app_base import HeliumAppBase
from heliumtools.misc.gather_data import apply_ROI

# Configuration des logs quand le programme est appelé localment
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class FiguresClass:
    """This class contains the figure that is shown to the user, on the right of the figure."""

    def __init__(self, fig, canvas):
        self.figures = fig
        self.canvas = canvas

    def update_plot(self, model, tab_number):
        self.Tab1Fig(model, tab_number)

    def Tab1Fig(self, model, i):
        """Define here the figure of your first tab (i= 0 a priori)."""
        self.figures[i].clf()
        ax = self.figures[i].add_subplot(111)
        x = np.linspace(min(model.x_range), max(model.x_range), 500)
        y = gaussian(x, model.mean, model.amplitude, model.std)
        ax.plot(x, y)
        ax.set_xlabel("My x axis")
        ax.set_ylabel(r"My Y axis  $\int \partial P dt$")
        ax.grid(True)
        ax.set_title(
            r"Gaussian function with width {:.2f}".format(model.std),
            fontsize="medium",
        )
        self.figures[i].tight_layout()
        self.canvas[i].draw()


class ExampleApplication(HeliumAppBase):
    def __init__(self, data=[], metadata={}, parameters={}):
        """
        Initialisation of the GUIs. Before initialisation of the GUI, we must create
        """
        # Define the model
        model = Model(data=data, metadata=metadata, parameters=parameters)
        tab_names = ["Gaussian Curve", "Sinc Curve", "Test"]
        super().__init__(model=model, tab_names=tab_names)
        # Connect the FigureClass to your widget
        self.figures_class = FiguresClass(
            self.tabs_fig_widget.figures, self.tabs_fig_widget.canvas
        )

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
        # Built the windows Main Layout from left and right parts.
        self._main_widget = QSplitter()
        self._main_widget.addWidget(self.parameter_widget_wrapper)
        self._main_widget.addWidget(self.right_widget)
        self.setCentralWidget(self._main_widget)
        # self.main_layout = QHBoxLayout()
        # self.main_layout.addLayout(self.layout_left)
        # self.main_layout.addLayout(self.layout_right)

        # self.setLayout(self.main_layout)

    def _setup_UI_left_part(self):
        """This method sets up the left part of the UI aka the parameter part."""
        self.left_widget = self.parameter_widget_wrapper

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
