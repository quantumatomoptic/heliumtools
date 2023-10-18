import matplotlib.pyplot as plt
import numpy as np
from math import ceil
import seaborn as sns


class FigureClass:
    """This class contains the figure that is shown to the user, on the right of the figure."""

    def __init__(self, fig, canvas):
        self.figure = fig
        self.canvas = canvas

    def update_plot(self, model):
        """Define here the figure."""
        self.figure.clf()
        self.ncols = 3
        self.nrows = ceil(len(model.all_delay) / self.ncols)
        for i, delay in enumerate(model.all_delay):
            ax = self.figure.add_subplot(self.nrows, self.ncols, i + 1)
            data = model.result[i]
            sns.heatmap(
                data.pivot(index="Vz1", columns="Vz2", values=model.value_to_plot),
                cmap=model.cmap,
                ax=ax,
                vmin=model.fig_vmin,
                vmax=model.fig_vmax,
            )
            ax.invert_yaxis()
            ax.set_title("Delay : {} us".format(delay), fontsize="medium")
        title = "Vz, Vz' correlations with "
        title += model.get_properties_for_title()
        self.figure.suptitle(
            title,
            fontsize="medium",
        )
        self.figure.tight_layout()
        self.canvas.draw()


if __name__ == "__main__":
    ma_fonction()
