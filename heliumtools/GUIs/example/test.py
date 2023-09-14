import sys
from PyQt5.QtWidgets import QApplication, QMainWindow
from widget import MonWidget


class MaFenetre(QMainWindow):
    def __init__(self):
        super().__init__()

        self.initUI()

    def initUI(self):
        self.setWindowTitle("Application à l'écoute des boutons")
        self.setGeometry(100, 100, 400, 200)

        # Créez une instance de votre widget personnalisé
        self.mon_widget = MonWidget()

        # Connectez le signal personnalisé du widget à une fonction de réaction
        self.mon_widget.bouton_presse.connect(self.reagir_au_bouton_presse)

        # Ajoutez le widget à la fenêtre principale
        self.setCentralWidget(self.mon_widget)

    def reagir_au_bouton_presse(self):
        # Fonction de réaction lorsque le bouton est pressé
        print("Bouton dans le widget pressé")


def main():
    app = QApplication(sys.argv)
    fenetre = MaFenetre()
    fenetre.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
