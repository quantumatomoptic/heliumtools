import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QLabel

# Crée une application Qt
app = QApplication(sys.argv)

# Crée une fenêtre principale
window = QMainWindow()
window.setWindowTitle('Ma Application PyQt5')
window.setGeometry(100, 100, 400, 200)  # Définit la position et la taille de la fenêtre

# Ajoute un label à la fenêtre
label = QLabel('Bonjour, PyQt5 !', window)
label.setGeometry(50, 50, 300, 100)  # Définit la position et la taille du label dans la fenêtre

# Affiche la fenêtre
window.show()

# Exécute l'application
sys.exit(app.exec_())
