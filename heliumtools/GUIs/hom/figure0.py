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
        columns_to_plot = [
            "BEC Std Arrival Time",
            "Number of Atoms",
            "Number of Atoms in ROI",
            "BEC Arrival Time",
            "dT of atoms in ROI",
            "Number of Atoms in supplementary ROI{}",
        ]
        ncols = 2
        nrows = ceil(len(columns_to_plot) / ncols)
        for index, column in enumerate(columns_to_plot):
            ax = self.figure.add_subplot(ncols, nrows, index + 1)
            if column not in model._all_bec_arrival_times.columns:
                continue
            sns.lineplot(
                data=model._all_bec_arrival_times,
                x=model._all_bec_arrival_times.index,
                y=column,
                ax=ax,
                alpha=0.8,
                hue="Sequence",
                palette="Dark2",
            )
            # coloriage en gris des zones où on a bien des cycles
            ax.fill_between(
                model._all_bec_arrival_times.index,
                np.min(model.bec_arrival_times[column])
                * np.ones(len(model._all_bec_arrival_times)),
                np.max(model.bec_arrival_times[column])
                * np.ones(len(model._all_bec_arrival_times)),
                alpha=0.25,
                color="Grey",
                label="Selection",
            )
            ax.grid(True, alpha=0.6)
            if index == 0:
                ax.legend(fontsize="x-small", ncols=2)
            else:
                ax.get_legend().remove()
            if column in model.filters:
                ax.plot(
                    model.bec_arrival_times.index,
                    model.filters[column][0] * np.ones(len(model.bec_arrival_times)),
                    ls="--",
                    alpha=0.8,
                    color="darkred",
                    label="Filter",
                )
                ax.plot(
                    model.bec_arrival_times.index,
                    model.filters[column][1] * np.ones(len(model.bec_arrival_times)),
                    ls="--",
                    alpha=0.8,
                    color="darkred",
                )
                ecart = np.abs(model.filters[column][1] - model.filters[column][0])
                bottom = max(
                    [
                        model.filters[column][0] - 2 * ecart,
                        np.min(model._all_bec_arrival_times[column]),
                    ]
                )
                top = min(
                    [
                        model.filters[column][1] + 2 * ecart,
                        max(model._all_bec_arrival_times[column]),
                    ]
                )
                ax.set_ylim(bottom=bottom, top=top)
        self.figure.suptitle(
            r"Stability of the selection : "
            + r"$T_{{BEC}} =  {:.3f} \pm {:.3f} $ ms ".format(
                model.bec_arrival_times["BEC Arrival Time"].mean(),
                model.bec_arrival_times["BEC Arrival Time"].std(),
            )
            + r"$N_{{pairs}} =  {:.0f} \pm {:.0f}$".format(
                model.bec_arrival_times["Number of Atoms in ROI"].mean(),
                model.bec_arrival_times["Number of Atoms in ROI"].std(),
            )
            + r" keeping {}/{} cycles.".format(
                len(model.bec_arrival_times), len(model._all_bec_arrival_times)
            ),
            fontsize="medium",
        )

        self.figure.tight_layout()
        self.canvas.draw()


if __name__ == "__main__":
    ma_fonction()
