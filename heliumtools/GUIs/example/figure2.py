import matplotlib.pyplot as plt
import numpy as np


def function(x, mean, amplitude, largeur):
    return np.where(
        abs((x - mean) / largeur) < 1,
        amplitude * (1 - ((x - mean) / largeur) ** 2) ** 2,
        0,
    )


class FigureClass:
    """This class contains the figure that is shown to the user, on the right of the figure."""

    def __init__(self, fig, canvas):
        self.figure = fig
        self.canvas = canvas

    def update_plot(self, model):
        """Define here the figure."""
        self.figure.clf()
        ax = self.figure.add_subplot(111)
        x = np.linspace(min(model.x_range), max(model.x_range), 500)
        y = function(x, model.mean, model.amplitude, model.std)
        ax.plot(x, y)
        ax.set_xlabel("My x axis")
        ax.set_ylabel(r"My Y axis")
        ax.grid(True)
        ax.set_title(
            r"Thomas Fermi profile with amplitude {:.2f}".format(model.amplitude),
            fontsize="medium",
        )
        self.figure.tight_layout()
        self.canvas.draw()
