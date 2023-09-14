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


class HeliumAppBase(QMainWindow):
    def __init__(self):
        super().__init__()
