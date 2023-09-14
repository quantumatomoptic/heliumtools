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
import logging, os
from heliumtools.GUIs.base.parameter_widget_base import ParameterBaseWidget
from heliumtools.GUIs.base.figures_in_tab_widget_base import TabsFigureBase

logger = logging.getLogger(__name__)


class HeliumAppBase(QMainWindow):
    def __init__(self, model, tab_names):
        super().__init__()
        self.model = model
        self.tab_names = tab_names
        # -- Create parameters widget
        self.parameter_widget_wrapper = ParameterBaseWidget(model=self.model)
        # -- Create the widgets with tabs and figures in it :
        self.tabs_fig_widget = TabsFigureBase(tab_names=tab_names)

    def set_icon(self, icon="atom.png"):
        path_to_icon = os.path.join("icons", icon)
        if os.path.exists(path_to_icon):
            self.setWindowIcon(QIcon(path_to_icon))
            return
        logger.warning(
            f" The icon you require {icon} is not in the standard folder. Trying to set the default icon. "
        )
        default_icon = os.path.join("icons", "atom.png")
        if os.path.exists(default_icon):
            self.setWindowIcon(QIcon(default_icon))

    def update_plot(self, tab_number):
        self.figure_class.update_plot(self.model, tab_number)
