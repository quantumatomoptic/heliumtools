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
Decription of app_base.py 

This class is designd to be the master class of your App. It dialogs with the
model you defined, the later being a attribut of your app class and with the
FigureClasses you defined in your file figure0, figure1 etc....
"""


from PyQt5.QtWidgets import (
    QApplication,
    QWidget,
    QLabel,
    QLineEdit,
    QPushButton,
    QVBoxLayout,
    QHBoxLayout,
    QSizePolicy,
    QMainWindow,
    QScrollArea,
    QSplitter,
    QAction,
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QIcon
import logging, os, importlib, sys
from heliumtools.GUIs.base.parameter_widget_base import ParameterBaseWidget
from heliumtools.GUIs.base.figures_in_tab_widget_base import TabsFigureBase
from heliumtools.GUIs.base.default_figure import DefaultFigureClass

logger = logging.getLogger(__name__)
import traceback


class HeliumAppBase(QMainWindow):
    def __init__(self, model, tab_names, file):
        super().__init__()
        self.model = model
        self.tab_names = tab_names
        self._app_file = file
        self._app_dir = os.path.dirname(self._app_file)
        self._set_default_UI_parameters()
        # - Use this method to create the User interface
        self.initUI()
        self.setup_file_menu()

    def initUI(self):
        """
        Main function that creates the User Interface.
        """
        # -- Create left and right part of the UI
        self._setup_UI_left_part()
        self._setup_UI_right_part()
        # -- Create the Central Widget of the App : a splitter containing
        #    left and right widgets.
        self._main_widget = QSplitter()
        self._main_widget.addWidget(self.left_widget)
        self._main_widget.addWidget(self.right_widget)
        self.setCentralWidget(self._main_widget)

    ##########################################################################
    ## ------ METHODS THAT ARE WRITTEN TO BE USED BY THE CHILD CLASS ------ ##
    ##########################################################################

    def add_user_widget(self, widget):
        """
        This function adds a widget to the left-down layout widget.
        """
        self.user_widgets_layout.addWidget(widget)

    def update_all_plots(self):
        """
        method that could be use by the child class to update all figures.
        """
        for i in range(len(self.tab_names)):
            self.update_plot(i)

    def update_plot(self, tab_number):
        """
        Method to be used by the user to update a the tab_number'th plot.
        """
        if type(tab_number) != int:
            logging.error(f"Tab number {tab_number} does not exist. ")
            return
        if (tab_number > len(self.figure_classes_imported)) or (
            tab_number < -len(self.figure_classes_imported)
        ):
            n = len(self.figure_classes_imported) - 1
            msg = f"The tab number should be an integer bewteen O and {n}"
            logging.error(msg)
            return
        try:
            self.figure_classes_imported[tab_number].update_plot(self.model)
        except Exception as e:
            msg = f"Updating plot number {tab_number} failed. Error type : {e}"
            logging.error(msg)
            traceback.print_exc()

    def update_model_with_user_parameters(self):
        """
        Method to be used by the user to update the model values according to
        the values entered by the user on the GUI.
        """
        self.model = self.parameter_widget.get_parameters_values(self.model)

    def set_icon(self, icon="atom.png"):
        """
        Method to be used by the user to set the app icon

        Parameters
        ----------
        icon : path or string, optional
            path to the app icon. The default is "atom.png".
        """
        if os.path.exists(icon):
            self.setWindowIcon(QIcon(icon))
            return
        path_to_icon = os.path.join("icons", icon)
        if os.path.exists(path_to_icon):
            self.setWindowIcon(QIcon(path_to_icon))
            return
        path_to_icon = os.path.join(os.path.dirname(__file__), "icons", icon)
        if os.path.exists(path_to_icon):
            self.setWindowIcon(QIcon(path_to_icon))
            return
        logger.warning(
            f" Did not found the {icon} you required. Trying to set the default icon. "
        )

    def save_parameters(self):
        self.model._save_parameters()

    ##########################################################################
    ## ---------- METHODS THAT ARE WRITTEN TO CONSTRUCT THE GUI ----------  ##
    ##########################################################################

    def _setup_UI_left_part(self):
        """
        This method sets up the left part of the UI. It is composed of two
        widgets in a vertical QSplitter(). The top widget is the parameter
        widget ParameterBaseWidget and the second one is a widget that contains
        a layout that should be filled with the user widgets.
        This method uses the setup_user_buttons() that should be defined in the
        child class.
        """
        # - define the left_widget : a splitter
        self.left_widget = QSplitter()
        self.left_widget.setOrientation(Qt.Vertical)

        # - Instanciate the parameter widget and add it to the Splitter
        self.parameter_widget = ParameterBaseWidget(model=self.model)
        self.left_widget.addWidget(self.parameter_widget)

        # - Create the widget that will contain the layout in which we will
        # add every user widgets
        self.user_widget_container = QWidget()
        self.user_widgets_layout = QVBoxLayout()
        # - Call the setup_user_buttons() method overiden by the child class to
        # create custom buttons
        self.setup_user_buttons()
        # Set layout to widget and you're done !
        self.user_widget_container.setLayout(self.user_widgets_layout)

        # Add the widget that contains user widgets to the UI
        self.left_widget.addWidget(self.user_widget_container)

    def _setup_UI_right_part(self):
        """
        This methods sets up the righ pannel of the GUI.
        This pannel is composed of only one (big) widgets defined in the
        TabsFigureBase class.
        Once created, we must connect the figure shown in the i-th tab to the
        figure that is actually defined by the user in the figurei.py script.
        We also need to map the 'Update Plot' button of each tab to the model.
        This is done connecting the signal sent by the TabsFigureBase
        (update_fig_signal, with the argument i ) to the update_plot function.
        """
        # -- Create the widgets with tabs and figures in it :
        self.tabs_fig_widget = TabsFigureBase(tab_names=self.tab_names)
        # -- Import figures
        self.import_and_connect_user_figures()
        # Now we must connect the signal that is emmitted each time
        # the button 'Update figure' is pressed in a tab
        self.tabs_fig_widget.update_fig_signal.connect(self.update_plot)
        # Define those widgets as the right part of the windows
        self.right_widget = self.tabs_fig_widget

    def setup_file_menu(self):
        file_menu = self.menuBar().addMenu("File")
        # Créez une action "Save"
        save_action = QAction("Save Parameters", self)
        save_action.triggered.connect(self.save_parameters)
        file_menu.addAction(save_action)

    def import_and_connect_user_figures(self):
        """
        Method that snoops into the folder of the app to find module named
        figureX.py, X ranging from 0 to the number of tab.
        If the methode find such a file, it instanciate the FigureCLass in it
        and adds it to the figure_classes_imported list. Each element of this
        list being the figure shown in the corresponding tab of the
        tabs_fig_widget widget.
        """
        # liste des FigureClass importées des fichiers figureX.py
        self.figure_classes_imported = []
        # on ajoute le chemin relatif du fichier (en fait pas nécessaire
        # mais je le laisse là au cas où)
        # sys.path.append(self._app_dir)
        for i in range(len(self.tab_names)):
            # -- import each class FigureClass from figureX.py
            try:
                module_name = f"figure{i}"
                mon_module = importlib.import_module(module_name)
                FigureClass = getattr(mon_module, "FigureClass")
                self.figure_classes_imported.append(
                    FigureClass(
                        self.tabs_fig_widget.figures[i], self.tabs_fig_widget.canvas[i]
                    )
                )
            except:
                # if we fail to import the module, we show the default figure.
                msg = f"Failed to load the FigureClass from the figure{i}.py"
                msg += f" file in {os.path.dirname(self._app_file)}. Please "
                msg += f"make sure you defined it well."
                logging.error(msg)
                self.figure_classes_imported.append(
                    DefaultFigureClass(
                        self.tabs_fig_widget.figures[i], self.tabs_fig_widget.canvas[i]
                    )
                )
            # Update the figure
            self.update_plot(i)
            # self.figure_classes_imported[i].update_plot(
            #     self.model,
            # )

    def _set_default_UI_parameters(self):
        self.default_icon = os.path.join(os.path.dirname(__file__), "icons", "atom.png")
        self.setWindowIcon(QIcon(self.default_icon))
        self.setWindowTitle("Example of an Application.")
        self.setGeometry(100, 100, 1000, 400)

    ##########################################################################
    ## -------------------- A PRIORI USELESS FUNCTIONS -------------------- ##
    ##########################################################################

    ### A priori useless functions that are here only so the App is self contained.
    def setup_user_buttons(self):
        """
        This function is supposed to be overide by the child class.
        Customize this method to create buttons related to custom callback functions.
        To do so :
            1. Create your button (or widget),
            2. Link it to a callback function that calls a function from your model
            3. Add this button to the 'button pannel' using the add_user_button(button) method
        """
        self.setup_user_buttons_default()

    def setup_user_buttons_default(self):
        """
        Default function that creates 2 useless buttons. Should be useless.
        """
        button1 = QPushButton("Useless button 1")
        button1.clicked.connect(self._no_callback_message)
        self.add_user_widget(button1)
        button2 = QPushButton("Press me button 2")
        button2.clicked.connect(self._no_callback_message)
        self.add_user_widget(button2)

    def _no_callback_message(self):
        """
        Warning message is shown when this message is called

        """
        msg = "The button you pressed has no callback. To create custom callbacks, please :\n"
        msg += "   1. Create your button (or widget), for example using 'button=QPushButton('my button')' \n"
        msg += "   2. link this button to a callback function (being for example: 'my_callback')."
        msg += "For example, you could use 'button.clicked.connect(my_callback). In your callback function, backreacts using methods of your model or updating a figure app using the update_plot(ith fig) defined in HeliumAppBase.\n"
        msg += "   3. Add this button to the 'button pannel' using the add_user_button(button) method."
        logging.warn(msg)
