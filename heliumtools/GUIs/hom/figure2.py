import matplotlib.pyplot as plt
import numpy as np
from math import ceil, sqrt
import seaborn as sns
import pandas as pd
from heliumtools.misc.gather_data import apply_ROI


class FigureClass:
    """This class contains the figure that is shown to the user, on the right of the figure."""

    def __init__(self, fig, canvas):
        self.figure = fig
        self.canvas = canvas

    def update_plot(self, model):
        """Define here the figure."""
        self.figure.clf()
        df = pd.concat(model.result)
        # on enlève les points sans Bragg
        df = df[df["Beam splitter delay (us)"] > 0]
        # on ne garde que les classes séparées de vBragg
        df = df[df["Vz1-Vz2"] == model.bragg_speed]
        # puis on ne garde que les ourbes pas trop loin du 'bon' référentiel
        df = apply_ROI(
            df,
            {
                "Vz1+Vz2": [
                    0 - model.oneD_frame_gap,
                    0 + model.oneD_frame_gap,
                ],
            },
        )
        # liste des courbes
        plot_by = "Vz1+Vz2"
        list_of_curves = df[plot_by].unique()
        n_plots = ceil(len(list_of_curves) / model.oneD_curves_per_plot)
        self.ncols = round(sqrt(n_plots))
        self.nrows = ceil(n_plots / self.ncols)
        newarr = np.array_split(np.array(list_of_curves), n_plots)
        for i, vz in enumerate(newarr):
            ax = self.figure.add_subplot(self.nrows, self.ncols, i + 1)
            data = df[df[plot_by].isin(vz)]
            sns.lineplot(
                data=data,
                x="Beam splitter delay (us)",
                y=model.value_to_plot,
                palette="Dark2",
                hue=plot_by,
                style=plot_by,
                ax=ax,
            )
            ax.legend(fontsize="x-small", ncols=2)
        title = "Looking for HOM dip with "
        title += model.get_properties_for_title()
        title += " $v_{{Bragg}} = {}$".format(model.bragg_speed)
        self.figure.suptitle(
            title,
            fontsize="medium",
        )
        self.figure.tight_layout()
        self.canvas.draw()


if __name__ == "__main__":
    ma_fonction()
