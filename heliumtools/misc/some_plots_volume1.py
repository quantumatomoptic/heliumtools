#!/usr/bin/env python
# -*- mode:Python; coding: utf-8 -*-

# ----------------------------------
#  Created on the Tue Jun 28 2022
#  Author : Victor
#  Copyright (c) 2022 - Helium1@LCF
# ----------------------------------
#
"""
Content of some_plots_volume1

Quelques plots pour laisser de l'espace au script principal.
"""

import matplotlib.pyplot as plt
from math import ceil
import seaborn as sns
import numpy as np


def histogram_plots_with_inset(df, X1, X2, Y, X1_label="", X2_label="", title=""):
    """_summary_

    Parameters
    ----------
    df : pandas dataframe
        avec les données à plotter
    X1 : str
        _description_
    X2 : str
        _description_
    Y : str
        _description_
    """
    if X1_label == "":
        X1_label = X1
    if X2_label == "":
        X2_label = X2
    fig, ax1 = plt.subplots()
    sigma0 = df[X1].std()
    mean0 = df[X1].mean()
    plt.hist(
        df[X1],
        bins=100,
        color="steelblue",
        label="Statistical properties\n mean = {:.3f} \n σ ={:.3f} ms ".format(
            mean0, sigma0
        ),
    )
    plt.legend(loc="upper right")
    plt.xlabel(X1_label)
    plt.ylabel(Y)
    if title != "":
        plt.title(title)
    # Add the number of atoms
    [[a1, b1], [a2, b2]] = ax1.get_position().get_points()
    height = 0.2
    width = 0.3
    alpha = 0.5
    color = "olivedrab"
    ax2 = fig.add_axes([a1 + 0.05, b2 - height - 0.05, width, height])
    ax2.hist(df[X2], color=color, alpha=alpha, bins=30)
    ax2.tick_params(labelsize=8)
    ax2.patch.set_alpha(alpha)
    plt.title(f"{X2} hist", fontsize=9)
    plt.show()


def stability_of_sequence(
    bec_arrival_times,
    columns_to_plot=[
        "Number of Atoms in ROIfit",
        "Number of Atoms in ROI",
        "BEC Arrival Time",
        "BEC Std Arrival Time",
    ],
    selec_bec_arrival_times=None,
    filters={},
    **kwargs,
):
    """Function that shows the stability of the sequence in each column of the bec arrival_times dataframe. The dashed red light show the filters used for the selected cycles and the grey area shows the global region of the made selection.

    Parameters
    ----------
    bec_arrival_times : df
        pandas dataframe with arrival times
    columns_to_plot : list, optional
        columns that on wants to check, by default [ "Number of Atoms in ROIfit", "Number of Atoms in ROI", "BEC Arrival Time", "BEC Std Arrival Time", ]
    selec_bec_arrival_times : pandas dataframe, optional
        pndas dataframe with cycles chosen using filters, by default None
    filters : dict, optional
        filters use to select the non pathological filters. Elements should be tuples or lists with min and max, by default {}
    """
    ncols = 2
    nrows = ceil(len(columns_to_plot) / ncols)
    fig, axes = plt.subplots(ncols=ncols, figsize=(11, 3 * nrows), nrows=nrows)
    if selec_bec_arrival_times is None:
        selec_bec_arrival_times = bec_arrival_times
    for i, column in enumerate(columns_to_plot):
        ax = axes.flatten()[i]
        sns.lineplot(
            data=bec_arrival_times,
            x=bec_arrival_times.index,
            y=column,
            ax=ax,
            hue="Sequence",
            palette="Dark2",
        )
        ax.fill_between(
            bec_arrival_times.index,
            np.min(selec_bec_arrival_times[column]) * np.ones(len(bec_arrival_times)),
            np.max(selec_bec_arrival_times[column]) * np.ones(len(bec_arrival_times)),
            alpha=0.25,
            color="Grey",
            label="Selection",
        )

        if column in filters:
            ax.plot(
                bec_arrival_times.index,
                filters[column][0] * np.ones(len(bec_arrival_times)),
                ls="--",
                alpha=0.8,
                color="darkred",
                label="Filter",
            )
            ax.plot(
                bec_arrival_times.index,
                filters[column][1] * np.ones(len(bec_arrival_times)),
                ls="--",
                alpha=0.8,
                color="darkred",
            )
            ecart = np.abs(filters[column][1] - filters[column][0])
            bottom = max(
                [filters[column][0] - 2 * ecart, np.min(bec_arrival_times[column])]
            )
            top = min([filters[column][1] + 2 * ecart, max(bec_arrival_times[column])])
            ax.set_ylim(bottom=bottom, top=top)

    fig.suptitle(
        r"Stability of the selection : "
        + r"$T_{{BEC}} =  {:.3f} \pm {:.3f} $ ms ".format(
            selec_bec_arrival_times["BEC Arrival Time"].mean(),
            selec_bec_arrival_times["BEC Arrival Time"].std(),
        )
        + r"$N_{{pairs}} =  {:.0f} \pm {:.0f}$".format(
            selec_bec_arrival_times["Number of Atoms in ROI"].mean(),
            selec_bec_arrival_times["Number of Atoms in ROI"].std(),
        )
        + r" keeping {}/{} cycles.".format(
            len(selec_bec_arrival_times), len(bec_arrival_times)
        )
    )

    plt.tight_layout()
    plt.show()
