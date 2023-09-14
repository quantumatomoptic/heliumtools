from PyQt5.QtWidgets import (
    QWidget,
    QPushButton,
    QVBoxLayout,
    QHBoxLayout,
    QSizePolicy,
    QTabWidget,
)
import logging
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import numpy as np
from PyQt5.QtCore import pyqtSignal


class TabsFigureBase(QTabWidget):
    update_fig_signal = pyqtSignal(int)

    def __init__(self, tab_names=["Gaussian Curve"]):
        super().__init__()
        self.tab_names = tab_names
        self.set_up_tabs()

    def set_up_tabs(self):
        self.tab_widget_list = []
        self.figures = []
        self.canvas = []
        for i, tab_name in enumerate(self.tab_names):
            onglet_widget, figure, canvas = self.create_tab(tab_number=i)
            self.canvas.append(canvas)
            self.figures.append(figure)
            self.tab_widget_list.append(onglet_widget)
            self.addTab(onglet_widget, tab_name)  # base function of a QTabWidget

    def create_tab(self, tab_number):
        # Créez un widget pour l'onglet
        onglet_widget = QWidget()
        # on va faire un layout pour mettre :
        layout = QVBoxLayout()
        # --1 la figure
        figure = Figure()
        canvas = FigureCanvas(figure)
        ax = figure.add_subplot(111)
        x = np.linspace(-1, 1, 200)
        ax.plot(x, x ** (tab_number + 1))
        canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        # qu'on ajoute au layout widget de l'onglet
        layout.addWidget(canvas)
        # --2 push button 'Update Figure'
        button = QPushButton("Update Figure")
        button.clicked.connect(lambda: self.signal_to_update_fig_callback(tab_number))
        layout.addWidget(button)
        # --3 push button 'Open in matplotlib'
        button = QPushButton("Open in Matplotlib")
        button.clicked.connect(
            lambda: self.open_in_matplotlib_window_callback(tab_number)
        )
        layout.addWidget(button)
        # Set layout to widget and you're done !
        onglet_widget.setLayout(layout)

        return onglet_widget, figure, canvas

    def open_in_matplotlib_window_callback(self, tab_number=0):
        """Callback function. When the 'Open in Matplotlib' button of the i-th tab is pressed, we open it in a matplotlib figure."""
        dummy = plt.figure()
        new_manager = dummy.canvas.manager
        new_manager.canvas.figure = self.figures[tab_number]
        self.figures[tab_number].set_canvas(new_manager.canvas)
        plt.show()

    def signal_to_update_fig_callback(self, tab_number):
        logging.debug(f"Button pushed on tab {tab_number}")
        self.update_fig_signal.emit(tab_number)
