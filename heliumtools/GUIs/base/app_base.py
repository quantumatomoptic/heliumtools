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
)
from PyQt5.QtGui import QIcon
import logging, os, importlib, sys
from heliumtools.GUIs.base.parameter_widget_base import ParameterBaseWidget
from heliumtools.GUIs.base.figures_in_tab_widget_base import TabsFigureBase
from heliumtools.GUIs.base.default_figure import DefaultFigureClass

logger = logging.getLogger(__name__)


class HeliumAppBase(QMainWindow):
    def __init__(self, model, tab_names, file):
        super().__init__()
        self.model = model
        self.tab_names = tab_names
        self._app_file = file
        self._app_dir = os.path.dirname(self._app_file)
        # -- Create parameters widget
        self.parameter_widget = ParameterBaseWidget(model=self.model)
        # -- Create the widgets with tabs and figures in it :
        self.tabs_fig_widget = TabsFigureBase(tab_names=tab_names)
        # -- Import figures
        self.import_and_connect_user_figures()

    def import_and_connect_user_figures(self):
        self.figure_models = []  # liste des FigureClass importées
        # on ajoute le chemin relatif du fichier (en fait pas nécessaire)
        # sys.path.append(self._app_dir)
        for i in range(len(self.tab_names)):
            # -- import each class
            try:
                module_name = f"figure{i}"
                mon_module = importlib.import_module(module_name)
                FigureClass = getattr(mon_module, "FigureClass")
                self.figure_models.append(
                    FigureClass(
                        self.tabs_fig_widget.figures[i], self.tabs_fig_widget.canvas[i]
                    )
                )
            except:
                logging.error(
                    f" Failed to load the FigureClass from the figure{i}.py file in  {self._app_file}. Please make sure you defined it well."
                )
                self.figure_models.append(
                    DefaultFigureClass(
                        self.tabs_fig_widget.figures[i], self.tabs_fig_widget.canvas[i]
                    )
                )
            self.figure_models[i].update_plot(
                self.model,
            )

    def set_icon(self, icon="atom.png"):
        path_to_icon = os.path.join("icons", icon)
        if os.path.exists(path_to_icon):
            self.setWindowIcon(QIcon(path_to_icon))
            return
        logger.warning(
            f" The icon you require {icon} is not in the standard folder. Trying to set the default icon. "
        )
        default_icon = os.path.join(os.path.dirname(__file__), "icons", "atom.png")
        if os.path.exists(default_icon):
            self.setWindowIcon(QIcon(default_icon))

    def update_plot(self, tab_number):
        self.figure_class.update_plot(self.model, tab_number)
