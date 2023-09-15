from PyQt5.QtWidgets import (
    QWidget,
    QLabel,
    QLineEdit,
    QVBoxLayout,
    QHBoxLayout,
    QSizePolicy,
    QScrollArea,
)
from PyQt5.QtCore import Qt

import snoop


class ParameterBaseWidget(QWidget):
    def __init__(self, model):
        super().__init__()
        self._built_parameter_widget(model)

    def _built_parameter_widget(self, model):
        parameter_widget_layout = QVBoxLayout(self)
        # title of the coluns
        label = QLabel("<html><b> Parameters </b></html>")
        label.setAlignment(Qt.AlignCenter)  # Alignement centré
        parameter_widget_layout.addWidget(label)
        # Define a scroll area for parameters
        self.scroll_area = QScrollArea()
        self.scroll_widget = QWidget()
        self.scroll_layout = QVBoxLayout()
        self.scroll_widget.setLayout(self.scroll_layout)
        # liste contenant les QLabel des paramètres
        self.labels_of_parameters = []
        # liste contenant les QLine edits des paramètres
        self.values_of_parameters = []
        # on parourt la liste des paramètres du modèle
        for name, value in model.get_parameters_str().items():
            self.labels_of_parameters.append(QLabel(name))
            self.values_of_parameters.append(QLineEdit())
            # set the default value in it
            self.values_of_parameters[-1].setText(str(value))
            # creat a horizontal layout foqvalue.text()r the name and the value of the parameter
            small_layout = QHBoxLayout()
            small_layout.addWidget(self.labels_of_parameters[-1])
            small_layout.addWidget(self.values_of_parameters[-1])
            # and store it into the vertical left layout
            self.scroll_layout.addLayout(small_layout)
        self.scroll_area.setWidget(self.scroll_widget)
        self.scroll_area.setWidgetResizable(True)
        parameter_widget_layout.addWidget(self.scroll_area)

    def get_parameters_values(self, model):
        for qlabel, qvalue in zip(self.labels_of_parameters, self.values_of_parameters):
            model_value = model.update_parameter(qlabel.text(), qvalue.text())
            qvalue.setText(model_value)
        return model
