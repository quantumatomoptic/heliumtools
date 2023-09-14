import matplotlib.pyplot as plt
import numpy as np


def gaussian(x, mean, amplitude, largeur):
    return amplitude * np.exp(-((x - mean) ** 2) / 2 / largeur**2)


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
        y = gaussian(x, model.mean, model.amplitude, model.std)
        ax.plot(x, y)
        ax.set_xlabel("My x axis")
        ax.set_ylabel(r"My Y axis  $\int \partial P dt$")
        ax.grid(True)
        ax.set_title(
            r"Gaussian function with width {:.2f}".format(model.std),
            fontsize="medium",
        )
        self.figure.tight_layout()
        self.canvas.draw()


if __name__ == "__main__":
    ma_fonction()
