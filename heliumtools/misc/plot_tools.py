#!/usr/bin/env python
# -*- mode:Python; coding: utf-8 -*-

# ----------------------------------
# Created on the 16-5-2023 by Victor
# Developped by Victor, ...
# Copyright (c) 2023 - Helium1@LCF
# ----------------------------------
#
"""
Content of plot_tools.py
-----------------------------

Please document your code ;-).
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap, LinearSegmentedColormap


def create_linear_colormap(
    RGB1=(1, 0, 0),
    RGB2=(0, 1, 0),
    pivot=True,
    n_points=256,
    RGB3=(0, 0, 1),
    n_points1=None,
    n_points2=None,
    name="MyWondefullColormap",
):
    """Renvoie une colormap matplotlib interpolée linéairement entre deux couleurs RGB1 et RGB2n avec éventuellement une troisième couleur pivot.

    Parameters
    ----------
    RGB1 : tuple, optional
        Couleur RGB pour les plus petites valeurs de la colormap, by default (1,0,0)
    RGB2 : tuple, optional
        Couleur RGB pour les plus grandes valeurs de la colormap, by default (0,1,0)
    pivot : bool, optional
        S'il y a une troisième valeur pour la colormap, by default True
    n_points : int, optional
        le nombre de points pour créer la colormap, by default 256
    RGB3 : tuple, optional
        la valeur RGB du troisième point (pivot, au milieu de la colormap) , by default (0,0,1)
    n_points1 : int, optional
        le nombre de points pour le premier segment, dans le cas où il y a une couleur pivot, by default None
    n_points2 : _type_, optional
        le nombre de points pour le second segment, dans le cas où il y a une couleur pivot, by default None

    Returns
    -------
    matplotlib.ListedColormap
        Colormap interpolée entre RGB1 et RGB2, voir RGB3 si pivot est vrai.
    """
    if pivot:
        if n_points1 == None:
            n_points1 = n_points
        list1 = np.array(
            [np.linspace(RGB1[i], RGB3[i], n_points1) for i in range(3)]
            + [np.ones(n_points1)]
        ).T
        if n_points2 == None:
            n_points2 = n_points
        list2 = np.array(
            [np.linspace(RGB3[i], RGB2[i], n_points2) for i in range(3)]
            + [np.ones(n_points2)]
        ).T
        newcolors = np.vstack((list1, list2))
    else:
        newcolors = np.array(
            [np.linspace(RGB1[i], RGB2[i], n_points) for i in range(3)]
            + [np.ones(n_points)]
        ).T
    return ListedColormap(newcolors, name=name)


blue_iogs_rgb = (10 / 255, 50 / 255, 80 / 255)
blue_iogs_cmjn = (0 / 255, 56 / 255, 101 / 255)
red_iogs_rgb = (255 / 255, 150 / 255, 1 / 255)
red_iogs_cmjn = (254 / 255, 80 / 255, 0)
white = (1, 1, 1)
black = (0, 0, 0)

cmap_iogs_cmjn = create_linear_colormap(
    RGB1=blue_iogs_cmjn,
    RGB2=white,
    RGB3=red_iogs_cmjn,
    n_points=300,
    pivot=True,
    name="IOGS_CMJN",
)

cmap_iogs_cmjn_white_pivot = create_linear_colormap(
    RGB1=blue_iogs_cmjn,
    RGB2=red_iogs_cmjn,
    RGB3=white,
    n_points=300,
    pivot=True,
    name="IOGS_CMJN_WHITE_PIVOT",
)
cmap_iogs_rgb = create_linear_colormap(
    RGB1=blue_iogs_rgb,
    RGB2=white,
    RGB3=red_iogs_rgb,
    n_points=300,
    pivot=True,
    name="IOGS_RGB",
)
