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
from heliumtools.misc.gather_data import apply_ROI
from scipy.optimize import curve_fit


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
    selec_bec_arrival_times = apply_ROI(bec_arrival_times, filters)
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


def show_2D_density_XY_with_fit(corr):
    peak1 = corr.atoms[(corr.atoms["Vz"] < -10) & (corr.atoms["Vz"] > -40)]
    peak2 = corr.atoms[(corr.atoms["Vz"] < 40) & (corr.atoms["Vz"] > 10)]
    from heliumtools.fit.gauss2D import Gauss2DFit, Gauss2D
    from heliumtools.fit.gauss1D import Gauss1DFit
    from heliumtools.fit.lorentz2D import Lorentz2DFit
    from heliumtools.fit.lorentz1D import Lorentz1DFit

    n_bins = 50

    fig, axes = plt.subplots(figsize=(12, 8), ncols=3, nrows=2)
    # Première colonne : densité
    ax = axes[0, 0]
    sns.histplot(
        data=peak1,
        x="Vx",
        y="Vy",
        ax=axes[0, 0],
        cmap="viridis",
        cbar=True,
        bins=(n_bins, n_bins),
    )
    hist_data, vx_data, vy_data = np.histogram2d(
        peak1["Vx"], peak1["Vy"], bins=n_bins, range=([-60, 18], [-60, 60])
    )
    vy_data = (vy_data[0:-1] + vy_data[1:]) / 2
    vx_data = (vx_data[0:-1] + vx_data[1:]) / 2
    VY, VX = np.meshgrid(vy_data, vx_data)
    g2Dfit = Gauss2DFit(x=(VX, VY), z=hist_data)
    g2Dfit.do_guess()
    g2Dfit.do_fit()
    fit_data = g2Dfit.eval(x=(VX, VY))
    l2Dfit = Lorentz2DFit(x=(VX, VY), z=hist_data)
    l2Dfit.do_guess()
    l2Dfit.do_fit()
    fit_data_lorentz = l2Dfit.eval(x=(VX, VY))
    popt1 = g2Dfit.popt
    axes[0, 0].set_title(
        "Peak 1 : -40 < Vz < -10 mm/s\n"
        + r"$\sigma_x=${:.0f} mm/s and $\sigma_y$={:.0f} mm/s".format(
            g2Dfit.popt[2], g2Dfit.popt[3]
        )
        + "\n"
        + "$c_x = ${:.1f} mm/s and $c_y$ ={:.1f} mm/s".format(
            g2Dfit.popt[4], g2Dfit.popt[5]
        )
    )
    ## Coupe 1D Selon X, Y est fixé au maximum de densité. On refait un fit mais 1D.
    minimumVY = np.argmin(np.abs(vy_data - g2Dfit.popt[5]))
    ax = axes[0, 1]
    g1Dfit = Gauss1DFit(x=vx_data, z=hist_data[:, minimumVY])
    g1Dfit.do_guess()
    g1Dfit.do_fit()
    l1Dfit = Lorentz1DFit(x=vx_data, z=hist_data[:, minimumVY])
    l1Dfit.do_guess()
    l1Dfit.do_fit()
    ax.plot(vx_data, hist_data[:, minimumVY], "o", alpha=0.6, label="Exp")
    ax.plot(vx_data, fit_data[:, minimumVY], "--", alpha=1, label="2D Gauss")
    ax.plot(vx_data, g1Dfit.eval(x=vx_data), "-", alpha=1, label="1D Gauss")
    ax.plot(vx_data, fit_data_lorentz[:, minimumVY], "-.", alpha=1, label="2D Lorentz")
    ax.plot(vx_data, l1Dfit.eval(x=vx_data), "-", alpha=1, label="1D Lorentz")
    ax.set_ylabel("# of atoms (a.u.)")
    ax.set_xlabel("Vx (mm/s)")
    ax.legend()
    ax.set_title(
        "1D cut along X at maximum density."
        + "\n"
        + "Gauss : $\sigma_x=${:.0f} mm/s ; $c_x = ${:.1f} mm/s".format(
            g1Dfit.popt[2], g1Dfit.popt[3]
        )
        + "\n"
        + "Lorentz : $\Gamma_x=${:.0f} mm/s ; $c_x = ${:.1f} mm/s".format(
            l1Dfit.popt[2], l1Dfit.popt[3]
        )
    )
    ## Coupe 1D Selon Y, X est fixé au maximum de densité. On refait un fit mais 1D.
    minimumVX = np.argmin(np.abs(vx_data - g2Dfit.popt[4]))
    ax = axes[0, 2]
    g1Dfit = Gauss1DFit(x=vy_data, z=hist_data[minimumVX, :])
    g1Dfit.do_guess()
    g1Dfit.do_fit()
    l1Dfit = Lorentz1DFit(x=vy_data, z=hist_data[minimumVX, :])
    l1Dfit.do_guess()
    l1Dfit.do_fit()
    ax.set_title("1D cut along Y at maximum density")
    ax.plot(vy_data, hist_data[minimumVX, :], "o", alpha=0.6)
    ax.plot(vy_data, fit_data[minimumVX, :], "--", alpha=1, label="2D Gaussian")
    ax.plot(vy_data, g1Dfit.eval(x=vy_data), "-", alpha=1, label="1D Gauss")
    ax.plot(vy_data, fit_data_lorentz[minimumVY, :], "-.", alpha=1, label="2D Lorentz")
    ax.plot(vy_data, l1Dfit.eval(x=vy_data), "-", alpha=1, label="1D Lorentz")
    ax.set_ylabel("# of atoms (a.u.)")
    ax.set_xlabel("Vy (mm/s)")
    ax.legend()
    ax.set_title(
        "1D cut along Y at maximum density."
        + "\n"
        + "Gauss : $\sigma_y=${:.0f} mm/s ; $c_y = ${:.1f} mm/s".format(
            g1Dfit.popt[2], g1Dfit.popt[3]
        )
        + "\n"
        + "Lorentz : $\Gamma_y=${:.0f} mm/s ; $c_y = ${:.1f} mm/s".format(
            l1Dfit.popt[2], l1Dfit.popt[3]
        )
    )

    sns.histplot(
        data=peak2,
        x="Vx",
        y="Vy",
        ax=axes[1, 0],
        cmap="viridis",
        cbar=True,
        bins=(n_bins, n_bins),
    )
    hist_data1, vx_data1, vy_data1 = np.histogram2d(
        peak2["Vx"], peak2["Vy"], bins=n_bins, range=([-60, 18], [-60, 60])
    )
    VY, VX = np.meshgrid(
        (vy_data1[0:-1] + vy_data1[1:]) / 2, (vx_data1[0:-1] + vx_data1[1:]) / 2
    )
    g2Dfit = Gauss2DFit(x=(VX, VY), z=hist_data1)
    g2Dfit.do_guess()
    g2Dfit.do_fit()
    popt2 = g2Dfit.popt
    axes[1, 0].set_title(
        "Peak 2 : 10 < Vz < 40 mm/s \n"
        + r"$\sigma_x=${:.0f} mm/s and $\sigma_y$={:.0f} mm/s".format(
            g2Dfit.popt[2], g2Dfit.popt[3]
        )
        + "\n"
        + "$c_x = ${:.1f} mm/s and $c_y$ ={:.1f} mm/s".format(
            g2Dfit.popt[4], g2Dfit.popt[5]
        )
    )
    ## Coupe 1D Selon X, Y est fixé au maximum de densité. On refait un fit mais 1D.
    minimumVY = np.argmin(np.abs(vy_data - g2Dfit.popt[5]))
    ax = axes[1, 1]
    g1Dfit = Gauss1DFit(x=vx_data, z=hist_data[:, minimumVY])
    g1Dfit.do_guess()
    g1Dfit.do_fit()
    l1Dfit = Lorentz1DFit(x=vx_data, z=hist_data[:, minimumVY])
    l1Dfit.do_guess()
    l1Dfit.do_fit()
    ax.plot(vx_data, hist_data[:, minimumVY], "o", alpha=0.6, label="Exp")
    ax.plot(vx_data, fit_data[:, minimumVY], "--", alpha=1, label="2D Gauss")
    ax.plot(vx_data, g1Dfit.eval(x=vx_data), "-", alpha=1, label="1D Gauss")
    ax.plot(vx_data, fit_data_lorentz[:, minimumVY], "-.", alpha=1, label="2D Lorentz")
    ax.plot(vx_data, l1Dfit.eval(x=vx_data), "-", alpha=1, label="1D Lorentz")
    ax.set_ylabel("# of atoms (a.u.)")
    ax.set_xlabel("Vx (mm/s)")
    ax.legend()
    ax.set_title(
        "1D cut along X at maximum density."
        + "\n"
        + "Gauss : $\sigma_x=${:.0f} mm/s ; $c_x = ${:.1f} mm/s".format(
            g1Dfit.popt[2], g1Dfit.popt[3]
        )
        + "\n"
        + "Lorentz : $\Gamma_x=${:.0f} mm/s ; $c_x = ${:.1f} mm/s".format(
            l1Dfit.popt[2], l1Dfit.popt[3]
        )
    )
    ## Coupe 1D Selon Y, X est fixé au maximum de densité. On refait un fit mais 1D.
    minimumVX = np.argmin(np.abs(vx_data - g2Dfit.popt[4]))
    ax = axes[1, 2]
    g1Dfit = Gauss1DFit(x=vy_data, z=hist_data[minimumVX, :])
    g1Dfit.do_guess()
    g1Dfit.do_fit()
    l1Dfit = Lorentz1DFit(x=vy_data, z=hist_data[minimumVX, :])
    l1Dfit.do_guess()
    l1Dfit.do_fit()
    ax.set_title("1D cut along Y at maximum density")
    ax.plot(vy_data, hist_data[minimumVX, :], "o", alpha=0.6)
    ax.plot(vy_data, fit_data[minimumVX, :], "--", alpha=1, label="2D Gaussian")
    ax.plot(vy_data, g1Dfit.eval(x=vy_data), "-", alpha=1, label="1D Gauss")
    ax.plot(vy_data, fit_data_lorentz[minimumVY, :], "-.", alpha=1, label="2D Lorentz")
    ax.plot(vy_data, l1Dfit.eval(x=vy_data), "-", alpha=1, label="1D Lorentz")
    ax.set_ylabel("# of atoms (a.u.)")
    ax.set_xlabel("Vy (mm/s)")
    ax.legend()
    ax.set_title(
        "1D cut along Y at maximum density."
        + "\n"
        + "Gauss : $\sigma_y=${:.0f} mm/s ; $c_y = ${:.1f} mm/s".format(
            g1Dfit.popt[2], g1Dfit.popt[3]
        )
        + "\n"
        + "Lorentz : $\Gamma_y=${:.0f} mm/s ; $c_y = ${:.1f} mm/s".format(
            l1Dfit.popt[2], l1Dfit.popt[3]
        )
    )

    plt.tight_layout()
    plt.show()


def oneD_density_fitted_plot(corr, miniVz=20, maxiVz=35, vperp_list=[10, 20, 30]):
    """Trace la densité 1D des paires selon z. 

    Parameters
    ----------
    corr : Correlation
    miniVz : float, optional
        vitesse minimal pour fit de la position de la paire, by default 20
    maxiVz : float, optional
        vitesse maximale pour le fit de la position de la paire, by default 35
    vperp_list : list, optional
        liste des vitesses transverses pour l'affichage, by default [10, 20, 30]

    Returns
    -------
    rien
    """
    peak1 = [-maxiVz, -miniVz]
    peak2 = [miniVz, maxiVz]
    boxZsize = corr.boxes["1"]["Vz"]["size"]

    def gaussian(x, mean, amplitude, standard_deviation):
        return amplitude * np.exp(-((x - mean) ** 2) / (2 * standard_deviation**2))

    #### Figure du nombre moyen d'atomes par boite :
    total_df = []
    corr.define_variable1(
        box="1", axe="Vz", type="position", name="Vz", min=-40, max=40, step=boxZsize
    )
    fig, ax = plt.subplots(figsize=(4.3, 3.9), dpi=110)
    for i, vperp_max in enumerate(vperp_list):
        my_box = {"Vperp": {"minimum": -1, "maximum": vperp_max}}
        atoms = corr.get_atoms_in_box(corr.atoms, my_box)
        hist, bins = np.histogram(atoms["Vz"], bins=np.arange(-40, 40, boxZsize))
        x = (bins[0:-1] + bins[1:]) / 2
        plt.scatter(
            x,
            hist / corr.n_cycles,
            alpha=0.3,
            label=r"$\Delta V_z={:.1f}, \, \Delta V_\perp={:.0f}$ mm/s ".format(
                boxZsize, vperp_max
            ),
        )

        ###### Fits de la paire 1
        try:
            hist1, bins1 = np.histogram(
                atoms["Vz"], bins=np.arange(peak1[0], peak1[1], boxZsize)
            )
            hist1 = hist1 / corr.n_cycles
            bin_centers1 = np.array(bins1[:-1] + np.diff(bins1) / 2)
            max_index = np.argmax(hist1)
            p0 = [
                bin_centers1[max_index],
                np.max(hist1),
                np.mean(hist1 * (bin_centers1 - bin_centers1[max_index])),
            ]
            popt, pcov_paire_1 = curve_fit(gaussian, bin_centers1, hist1, p0=p0)
            fit_absc = np.linspace(np.min(bin_centers1), np.max(bin_centers1), 50)
            plt.plot(fit_absc, gaussian(fit_absc, *popt), "--", color="C" + str(i))
            ax.text(
                0.8 * popt[0],
                popt[1],
                "{:.1f}".format(popt[0]),
                color="C" + str(i),
                ha="left",
                bbox=dict(facecolor="white", alpha=0.3, boxstyle="round"),
            )
        except Exception as exc:
            print(f"failed to fit pair 1. Error is {exc}")

        ###### Fits de la paire 2
        try:
            hist1, bins1 = np.histogram(
                atoms["Vz"], bins=np.arange(peak2[0], peak2[1], boxZsize)
            )
            hist1 = hist1 / corr.n_cycles
            bin_centers1 = np.array(bins1[:-1] + np.diff(bins1) / 2)
            max_index = np.argmax(hist1)
            # p0 = mean, amplitude, standard_deviation

            p0 = [
                bin_centers1[max_index],
                np.max(hist1),
                np.mean(hist1 * (bin_centers1 - bin_centers1[max_index]) ** 2),
            ]
            popt, pcov_paire_1 = curve_fit(
                gaussian, bin_centers1, hist1, p0=p0, maxfev=5000
            )
            fit_absc = np.linspace(np.min(bin_centers1), np.max(bin_centers1), 50)
            plt.plot(fit_absc, gaussian(fit_absc, *popt), "--", color="C" + str(i))
            ax.text(
                1.2 * popt[0],
                popt[1],
                "{:.1f}".format(popt[0]),
                color="C" + str(i),
                ha="left",
                bbox=dict(facecolor="white", alpha=0.3, boxstyle="round"),
            )
        except:
            print("failed to fit pair 2")

    fig.suptitle(
        r"$T_{{BEC}} = {:.3f}$ ms and $\vec{{V}}$=({}, {}, {}) mm/s".format(
            corr.bec_arrival_time,
            corr.ref_frame_speed["Vx"],
            corr.ref_frame_speed["Vy"],
            corr.ref_frame_speed["Vz"],
        ),
        fontsize="medium",
    )

    plt.ylabel("Mean number of atoms")
    plt.grid(True)
    plt.xlabel("Velocity along z (mm/s)")
    plt.legend(
        framealpha=0.4,
        fontsize="7",
    )
    # plt.savefig("densite_selon_z.png")
    plt.show()
