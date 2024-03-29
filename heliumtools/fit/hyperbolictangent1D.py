#!/usr/bin/env python3
# -*- mode: Python; coding: utf-8 -*-

"""
Author   : Charlie
Created  : 2021-05-03 10:49:51

Comments : implements a 1D hypetrangent fit
"""
# %% IMPORTS

# -- global
import numpy as np

# -- local
from heliumtools.fit.abstract import Abstract1DFit


# %% FUNCTIONS


def Hyperbolictangent1D(x, *p):
    return p[0] + p[1] * np.tanh((x - p[2]) / p[3])


# %% CLASS DEFINITION
class Hyperbolictangent1DFit(Abstract1DFit):
    """
    Class to fit datas : Hyperbolic Tangent

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
        self.name = "hyperbolic tangent"
        self.short_name = self.name
        self.category = "math"
        self.formula_help = "f(x) = p[0]+ p[1] * tanh(( x - p[2] ) / p[3])"
        self.parameters_help = "p = [offset, amplitude, center, slope]"
        self._version = "1.0"

    def _fitfunc(self, x, *p):
        return Hyperbolictangent1D(x, *p)

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
        zmin = np.min(z)
        zmax = np.max(z)
        xmin = np.min(x)
        xmax = np.max(x)
        A = zmax - zmin

        offset_guess = (zmax + zmin) / 2
        amp_guess = A / 2
        center_guess = (xmax + xmin) / 2
        slope_guess = -(xmax - xmin)
        # -- adapt to the current fit function
        p0 = [offset_guess, amp_guess, center_guess, slope_guess]

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
        offset, amplitude, center, slope = self.popt
        (
            offset_err,
            amplitude_err,
            center_err,
            slope_err,
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

        param = {
            "name": "center",
            "value": center,
            "error": center_err,
            "display": "%.3g",
            "unit": self.x_unit,
            "comment": "center",
        }
        values.append(param)

        param = {
            "name": "slope",
            "value": slope,
            "error": slope_err,
            "display": "%.3g",
            "unit": self.x_unit,
            "comment": "slope",
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
