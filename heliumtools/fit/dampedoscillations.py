#!/usr/bin/env python3
# -*- mode: Python; coding: utf-8 -*-
"""
Author   : Charlie
Created  : 2021-05-03 10:49:51

Comments : implements a 1D exponentially damped oscillations fit
"""
# %% IMPORTS

# -- global
import numpy as np

# -- local
from heliumtools.fit.abstract import Abstract1DFit


# %% FUNCTIONS


def DampedOscillation1D(x, *p):
    return p[0] + p[1] * np.sin(2 * np.pi * p[2] * x + p[3]) * np.exp(-x / p[4])


# %% CLASS DEFINITION
class DampedOscillation1DFit(Abstract1DFit):
    """
    Class to fit datas

    ----------
    USAGE :
    oneDfit = Hyperbolictangent1DFit(x=x, z=z)
    oneDfit.do_guess()
    oneDfit.guess = [*args] # to give you guess
    oneDfit.do_fit() #Now we do the fit using scipy.optimize.curve_fit

    ----------
    ATTRIBUTS
    popt = oneDfit.popt #optimal parameters
    perr = oneDfit.popt #errors on parameters
    pcov = oneDfit.pcov #covariance matrix
    formula_help = string, function formula
    parameters_help = list of string, help fo popt.

    ----------
    METHODS
    self.eval(x, params = None) : eval fit formula at coordinates 'x', with given set of 'params'.
    If params is not given, it uses the p0 found before.
    self.plot_fit_result() : plots the fit result, for a rapid check
    self.export_dic() : exports fit info and results as a python dictionnary
    self.export_json_str() : exports fit info and results as a json string
    self.do_guess() : guess optimal parameters

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # -- attributes specific to fit
        self.name = "damped oscillation"
        self.short_name = self.name
        self.category = "math"
        self.formula_help = "f(x) = p[0] "
        self.formula_help += "+ p[1] *sin(2*pi*p[2]*x+p[3]) * exp(-x/p[4])"
        self.parameters_help = (
            "p = [offset, amplitude, \n frequency, phase, damping rate]"
        )
        self._version = "1.0"

    def _fitfunc(self, x, *p):
        return DampedOscillation1D(x, *p)

    def do_guess(self):
        """guess fit parameters"""

        # -- check that the data and coordinates were provided
        if len(self.z) == 0 or len(self.x) == 0:
            return

        # -- get data
        z = self.z
        x = self.x

        # should be arrays
        x = np.asarray(x)
        z = np.asarray(z)

        # -- guess amplitude / offset / center / size

        # min, max, amp
        threshold = 0.05
        zmin = np.min(z)
        zmax = np.max(z)
        xmin = np.min(x)
        xmax = np.max(x)
        A = zmax - zmin

        offset_guess = (zmax + zmin) / 2
        amp_guess = A
        freq_guess = 3 / (xmax - xmin)
        freq_guess = 3 / (xmax - xmin)
        if len(x) > 3:
            # check if data are well separated. If not, keep 3.
            deltaX = x[1] - x[0]
            vector = x[1:-1] - x[0:-2] - deltaX
            # We check if datas are equally spaced, we can guess using fft
            if np.max(np.abs(vector)) < 0.0001:
                print("Entering")
                ft = np.fft.fft(z)
                N = len(z)
                T = N * deltaX
                argu = np.argmax(np.abs(ft[1 : int(N / 2)]))+1
                freq_guess = argu / T
        phase_guess = 0
        tau_guess = xmax - xmin
        # -- adapt to the current fit function
        p0 = [offset_guess, amp_guess, freq_guess, phase_guess, tau_guess]

        # save guess
        self.guess = p0

    def compute_values(self):
        """compute some physical values from the fit optimal parameters"""

        # -- check that the data and coordinates were provided
        if len(self.z) * len(self.x) * len(self.popt) == 0:
            return

        # -- get data
        z = self.z
        x = self.x
        zfit = self._fitfunc(x, *self.popt)

        # -- get fit results
        offset, amplitude, frequency, phase, damping_rate = self.popt
        (
            offset_err,
            amplitude_err,
            frequency_err,
            phase_err,
            damping_rate_err,
        ) = self.perr

        # -- init values list
        values = []

        # -- compute values
        # amplitude
        param = {
            "name": "amplitude",
            "value": amplitude,
            "error": amplitude_err,
            "display": "%.3g",
            "unit": self.z_unit,
            "comment": "amplitude",
        }
        values.append(param)

        # offset
        param = {
            "name": "offset",
            "value": offset,
            "error": offset_err,
            "display": "%.3g",
            "unit": self.z_unit,
            "comment": "offset",
        }
        values.append(param)

        # center

        if self.x_unit == "ms":
            frequency_unit = "kHz"
        elif self.x_unit == "us":
            frequency_unit = "MHz"
        elif self.x_unit == "s":
            frequency_unit = "Hz"
        else:
            frequency_unit = self.x_unit + "^-1"

        param = {
            "name": "frequency",
            "value": frequency,
            "error": frequency_err,
            "display": "%.3g",
            "unit": frequency_unit,
            "comment": "frequency",
        }
        values.append(param)
        # size
        param = {
            "name": "periode",
            "value": 1 / frequency,
            "error": frequency_err / (frequency ** 2),
            "display": "%.3g",
            "unit": self.x_unit,
            "comment": "periode",
        }
        values.append(param)
        # size
        param = {
            "name": "phase",
            "value": phase,
            "error": phase_err,
            "display": "%.3g",
            "unit": "rad",
            "comment": "phase at origin",
        }
        values.append(param)
        param = {
            "name": "damping_rate",
            "value": damping_rate,
            "error": damping_rate_err,
            "display": "%.3g",
            "unit": self.x_unit,
            "comment": "damping rate",
        }
        values.append(param)

        # -- spatial values in pixels
        # -- other
        # fit error
        fit_error = np.mean(np.sqrt((z - zfit) ** 2))
        param = {
            "name": "fit error",
            "value": fit_error,
            "display": "%.3g",
            "unit": "",
            "comment": "fit error = mean(sqrt((data - fit)**2)))",
        }
        values.append(param)

        # -- store
        self.values = values


# %%
