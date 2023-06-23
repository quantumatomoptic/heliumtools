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
