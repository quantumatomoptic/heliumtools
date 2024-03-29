#!/usr/bin/env python3
# -*- mode: Python; coding: utf-8 -*-

"""
Author   : Alexandre
Created  : 2021-05-03 10:49:51

Comments : implements a 1D Gauss fit
"""
# %% IMPORTS

# -- global
import numpy as np

# -- local
from heliumtools.fit.abstract import Abstract1DFit


# %% FUNCTIONS


def Gauss1DNoOffset(x, *p):
    return p[0] * np.exp(-((x - p[2]) ** 2) / 2 / p[1] ** 2)


# %% CLASS DEFINITION
class Gauss1DNoOffsetFit(Abstract1DFit):
    """
    Class to fit datas : 1D Gauss fit

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

        # -- attributes specific to 2D Gauss fit
        self.name = "gaussian"
        self.short_name = self.name
        self.category = "math"
        self.formula_help = "f(x) = 0 "
        self.formula_help += "+ p[0] * exp(-(x - p[2]) ** 2 / (2 * p[1]**2))"
        self.parameters_help = "p = [amplitude, size, center]"
        self._version = "1.0"

    def _fitfunc(self, x, *p):
        return Gauss1DNoOffset(x, *p)

    def do_guess(self):
        """guess fit parameters"""

        # -- check that the data and coordinates were provided
        if len(self.z) == 0 or len(self.x) == 0:
            return

        # -- get data
        z = self.z
        x = self.x  # this is a 2D fit !

        # should be arrays
        x = np.asarray(x)
        z = np.asarray(z)

        # -- guess amplitude / offset / center / size

        # min, max, amp
        threshold = 0.05
        zmin = np.min(z)
        zmax = np.max(z)
        A = zmax - zmin

        # offset, amplitude
        offset_guess = np.mean(z[z < zmin + threshold * A])
        amp_guess = zmax - offset_guess

        # filter to remove noise
        i_filter = (z - offset_guess) > threshold * amp_guess
        if len(i_filter) > 0:
            weights = z[i_filter] - offset_guess
            xf = x[i_filter]
            # center = center of mass
            c_guess = np.average(xf, weights=weights)
            # size = standard deviation
            s_guess = np.sqrt(np.average((xf - c_guess) ** 2, weights=weights))
        else:
            c_guess = np.mean(x)
            s_guess = 0.5 * (x.max() - x.min())

        # -- adapt to the current fit function
        p0 = [amp_guess, s_guess, c_guess]

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
        amplitude, size, center = self.popt
        amplitude_err, size_err, center_err = self.perr

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
            "comment": "gaussian amplitude",
        }
        values.append(param)

        # offset
        param = {
            "name": "offset",
            "value": 0,
            "error": 0,
            "display": "%.3g",
            "unit": self.z_unit,
            "comment": "gaussian offset",
        }
        values.append(param)

        # center
        param = {
            "name": "center",
            "value": center,
            "error": center_err,
            "display": "%.4f",
            "unit": self.x_unit,
            "comment": "gaussian center",
        }
        values.append(param)

        # size
        param = {
            "name": "size",
            "value": size,
            "error": size_err,
            "display": "%.3g",
            "unit": self.x_unit,
            "comment": "gaussian size",
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


# %% TESTS
if __name__ == "__main__":
    import matplotlib.pyplot as plt

    try:
        # -- generate noisy data
        x = np.linspace(-10, 10, 200)
        p = [0.5, 5.0, 0.5, 5]
        z = Gauss1D(x, *p)
        noise = (np.random.rand(*x.shape) - 0.5) * 2
        z += noise

        # -- Fit
        g1Dfit = Gauss1DFit(x=x, z=z)
        g1Dfit.do_guess()
        print(">> guess")
        print(g1Dfit.guess)
        g1Dfit.do_fit()
        xfit = np.linspace(x.min(), x.max(), 1000)
        zfit = g1Dfit.eval(xfit)

        print(">> popt")
        print(g1Dfit.popt)

        g1Dfit.compute_values()
        print(g1Dfit.export_json_str())

        plt.figure()
        plt.plot(x, z, "o")
        plt.plot(xfit, zfit)
        plt.show()
    except Exception as e:
        print(e)
