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
Once it is instanciate, the application gets back all parameters of the model from the function self.model.get_parameters(). 
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
from heliumtools.misc.parameter_model import GlobalModelForParameters
from heliumtools.misc.gather_data import apply_ROI

# Configuration des logs quand le programme est appelé localment
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def gaussian(x, mean, amplitude, standard_deviation):
    return amplitude * np.exp(-((x - mean) ** 2) / (2 * standard_deviation**2))


class Model(GlobalModelForParameters):
    """
    The model inherits the GlobalModelForParameters class from heliumtools.misc.
    """

    def __init__(self, parameters={}, data=[], metadata={}):
        self.mean = 2
        self.std = 4
        self.amplitude = 3
        self.x_range = [-20, 20]
        # The parameter dictionary given in the initialisation will overwrite
        # other parameters
        super().__init__(parameters=parameters)


class PlottingClasse:
    """This class contains the figure that is shown to the user, on the right of the figure."""

    def __init__(self, fig, ax, canvas):
        self.fig = fig
        self.ax = ax
        self.canvas = canvas

    def update_plot(self, model):
        self.ax.clear()
        x = np.linspace(min(model.x_range), max(model.x_range), 500)
        y = gaussian(x, model.mean, model.amplitude, model.std)
        self.ax.plot(x, y)
        self.ax.set_xlabel("My x axis")
        self.ax.set_ylabel(r"My Y axis  $\int \partial P dt$")
        self.ax.grid(True)
        self.ax.set_title(
            r"Gaussian function with width {:.2f}".format(model.std),
            fontsize="medium",
        )
        self.fig.tight_layout()
        self.canvas.draw()


class ExampleApplication(QWidget):
    def __init__(self, data=[], metadata={}, parameters={}):
        self.model = Model(data=data, metadata=metadata, parameters=parameters)
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.ax = self.figure.add_subplot(111)
        self.plot = PlottingClasse(self.figure, self.ax, self.canvas)
        super().__init__()
        self.initUI()

    def initUI(self):
        self.setWindowTitle("Density Visualizer.")
        self.setGeometry(100, 100, 1000, 400)
        try:
            self.setWindowIcon(
                QIcon(os.path.join("icons", "gauss.jpeg"))
            )  # Remplacez 'chemin_vers_votre_icone.ico' par le chemin de votre icône
        except:
            logger.warning(
                "WARNING : Our dear Gauss is missing. Please find him as fast as possible !!"
            )
        self._setup_UI_left_part()
        self._setup_UI_right_part()
        # Built the windows Main Layout from left and right parts.
        self.main_layout = QHBoxLayout()
        self.main_layout.addLayout(self.layout_left)
        self.main_layout.addLayout(self.layout_right)

        self.setLayout(self.main_layout)
        self.update_plot()

    def _setup_UI_left_part(self):
        """This method sets up the left part of the UI aka the parameter part."""
        # Left part of the windows : parameters
        self.layout_left = QVBoxLayout()
        # title of the coluns
        label = QLabel("<html><b> Parameters </b></html>")
        label.setAlignment(Qt.AlignCenter)  # Alignement centré
        self.layout_left.addWidget(label)
        # Define a scroll area for parameters
        self.scroll_area = QScrollArea()
        self.scroll_widget = QWidget()
        self.scroll_layout = QVBoxLayout()
        self.scroll_widget.setLayout(self.scroll_layout)
        self.labels_of_parameters = []  # liste contenant les QLabel des paramètres
        self.values_of_parameters = []  # liste contenant les QLine edits des paramètres
        # on parourt la liste des paramètres du modèle
        for name, value in self.model.get_parameters().items():
            self.labels_of_parameters.append(QLabel(name))
            self.values_of_parameters.append(QLineEdit())
            # set the default value in it
            self.values_of_parameters[-1].setText(str(value))
            # creat a horizontal layout foqvalue.text()r the name and the value of the parameter
            layout = QHBoxLayout()
            layout.addWidget(self.labels_of_parameters[-1])
            layout.addWidget(self.values_of_parameters[-1])
            # and store it into the vertical left layout
            self.scroll_layout.addLayout(layout)
        self.scroll_area.setWidget(self.scroll_widget)
        self.scroll_area.setWidgetResizable(True)
        self.layout_left.addWidget(self.scroll_area)

        ## now we add the button to update the figure
        self.update_figure_button = QPushButton("Update Figure")
        self.update_figure_button.clicked.connect(self.update_plot)  # connect action
        self.layout_left.addWidget(self.update_figure_button)

    def _setup_UI_right_part(self):
        """
        This function starts the right part of the GUI aka the figure part. It also instanciate the figure.
        If you modify the name of the PlottingClasse, do not forget to modify the
        """
        # Right part of the windows
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.layout_right = QVBoxLayout()
        self.layout_right.addWidget(self.canvas)

        ## now we add the button to update the figure
        self.open_in_matplotlib_window_button = QPushButton("Open in Matplotlib")
        self.open_in_matplotlib_window_button.clicked.connect(
            self.open_in_matplotlib_window
        )  # connect action
        self.layout_right.addWidget(self.open_in_matplotlib_window_button)
        # and plot now

    def update_plot(self):
        for qlabel, qvalue in zip(self.labels_of_parameters, self.values_of_parameters):
            model_value = self.model.update_parameter(qlabel.text(), qvalue.text())
            qvalue.setText(model_value)
        self.plot.update_plot(self.model)

    def open_in_matplotlib_window(self):
        """
        Callback function that is called when the user clicks on the 'Copy to Clip Board' button.
        """
        dummy = plt.figure()
        new_manager = dummy.canvas.manager
        new_manager.canvas.figure = self.figure
        self.figure.set_canvas(new_manager.canvas)
        plt.show()
        # # Sauvegarder la figure en tant qu'image PNG temporaire
        # temp_file = "temp_figure.png"
        # self.figure.savefig(temp_file, dpi=300, bbox_inches="tight")
        # # Charger l'image temporaire avec QImage
        # image = QImage(temp_file)
        # # Convertir QImage en QPixmap
        # pixmap = QPixmap.fromImage(image)
        # # Copier le QPixmap dans le presse-papiers
        # clipboard = QApplication.clipboard()
        # clipboard.setPixmap(pixmap)
        # # Supprimer le fichier temporaire
        # del image
        # os.remove(temp_file)


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
