from PyQt5.QtCore import QObject, pyqtSignal
from PyQt5.QtWidgets import QPushButton, QVBoxLayout, QWidget


class MonWidget(QWidget):
    # Créez un signal personnalisé pour le bouton pressé
    bouton_presse = pyqtSignal(int)

    def __init__(self):
        super().__init__()

        layout = QVBoxLayout()
        self.bouton = QPushButton("Bouton dans le widget")
        layout.addWidget(self.bouton)
        self.setLayout(layout)

        # Connectez le clic du bouton au signal personnalisé
        self.bouton.clicked.connect(self.emit_bouton_presse)

    def emit_bouton_presse(self):
        # Émettez le signal personnalisé lorsque le bouton est pressé
        self.bouton_presse.emit(2)
