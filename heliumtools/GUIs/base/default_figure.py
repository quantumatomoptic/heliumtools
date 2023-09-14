import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os


class DefaultFigureClass:
    """This class contains the figure that is shown to the user, on the right of the figure."""

    def __init__(self, fig, canvas):
        self.figure = fig
        self.canvas = canvas

    def update_plot(self, model):
        """Define here the figure."""
        self.figure.clf()
        file = os.path.join(os.path.dirname(__file__), "icons", "where.png")
        image = mpimg.imread(file)
        # Afficher l'image dans la figure
        ax = self.figure.add_subplot(111)
        ax.imshow(image)
        # Ajouter un titre à la figure (facultatif)
        ax.set_title(
            "One of your script is missing ! Check out the terminal.", fontsize="medium"
        )
        self.canvas.draw()


if __name__ == "__main__":
    ma_fonction()
