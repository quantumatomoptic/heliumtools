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
from heliumtools.misc.gather_data import (
    export_data_set_to_pickle,
    apply_ROI,
    load_atoms,
)
from heliumtools.correlations import Correlation

# Configuration des logs quand le programme est appelé localment
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class Model(GlobalModelForParameters):
    """
    The model inherits the GlobalModelForParameters class from heliumtools.misc.
    """

    def __init__(self, parameters={}, data=[], metadata={}):
        self.arrival_time = 307.5
        self.inertial_frame = [0, 0, 93]
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
        self.filters = {
            "BEC Arrival Time": [307.25, 307.75],
            "Number of Atoms in ROIfit": [800, 1500],
            "Number of Atoms": [3300, 5500],
            "Number of Atoms in ROI": [100, 280],
        }

        # The parameter dictionary given in the initialisation will overwrite
        # other parameters
        super().__init__(parameters=parameters)

    def compute_correlations(self):
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
            min=-self.var_max,
            max=-self.var_min,
            step=self.var_step,
        )
        self.corr.compute_correlations()
        logging.info("Computation done.")


class PlottingClasse:
    """This class contains the figure that is shown to the user, on the right of the figure."""

    def __init__(self, fig, ax, canvas):
        self.fig = fig
        self.ax = ax
        self.canvas = canvas

    def update_plot(self, model):
        self.ax.clear()
        self.fig.clf()
        self.ax = self.fig.add_subplot(111)

        sns.heatmap(
            model.corr.result.pivot(
                index="Vz1", columns="Vz2", values=model.value_to_plot
            ),
            cmap=model.cmap,
            ax=self.ax,
            vmin=model.fig_vmin,
            vmax=model.fig_vmax,
        )
        self.ax.invert_yaxis()
        title = (
            r"Correlations with $T_{{BEC}}$ = {} ms and $V_{{ref}}$={} mm/s. ".format(
                model.arrival_time, model.inertial_frame
            )
        )
        title += "\n"
        title += r"$\Delta V_z =$ {} mm/s ; $\Delta V_x $= {} mm/s,  $\Delta V_y $= {} mm/s".format(
            model.boxZsize, model.boxXsize, model.boxYsize
        )
        self.ax.set_title(
            title,
            fontsize="medium",
        )
        self.fig.tight_layout()
        self.canvas.draw()


class CorrelationApp(QWidget):
    def __init__(self, data=[], metadata={}, parameters={}):
        self.model = Model(data=data, metadata=metadata, parameters=parameters)
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.ax = self.figure.add_subplot(111)
        self.plot = PlottingClasse(self.figure, self.ax, self.canvas)
        super().__init__()
        self.initUI()

    def initUI(self):
        self.setWindowTitle("Correlation Visualizer.")
        self.setGeometry(100, 100, 1000, 700)
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
        self.compute_correlations_callback()

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

        ## -- -- -- Buttons
        # -- Update figure button
        self.compute_correlations_button = QPushButton("Compute correlations")
        self.compute_correlations_button.clicked.connect(
            self.compute_correlations_callback
        )  # connect action
        self.layout_left.addWidget(self.compute_correlations_button)

        # -- Update figure button
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

    def update_parameters(self):
        for qlabel, qvalue in zip(self.labels_of_parameters, self.values_of_parameters):
            model_value = self.model.update_parameter(qlabel.text(), qvalue.text())
            qvalue.setText(model_value)

    def compute_correlations_callback(self):
        self.update_parameters()
        self.model.compute_correlations()
        self.plot.update_plot(self.model)

    def update_plot(self):
        self.update_parameters()
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
    density_app = CorrelationApp(data=data, metadata=metadata, parameters=parameters)
    density_app.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    HOM_FOLDER = "/mnt/manip_E/2023/09/HOM/02"
    bec_arrival_times = pd.read_pickle(
        os.path.join(HOM_FOLDER, "bec_arr_time_BSdelay_0000_us.pkl")
    ).reset_index(drop=True)
    atoms = pd.read_pickle(
        os.path.join(HOM_FOLDER, "atoms_BSdelay_0000_us.pkl")
    ).reset_index(drop=True)
    parameters = {}
    main(atoms, bec_arrival_times, parameters)
