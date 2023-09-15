import matplotlib.pyplot as plt
import numpy as np


def function(x, mean, amplitude, largeur):
    return np.sinc((x - mean) / largeur) * amplitude


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
            r"Sinc function with width {:.2f}".format(model.std),
            fontsize="medium",
        )
        self.figure.tight_layout()
        self.canvas.draw()
