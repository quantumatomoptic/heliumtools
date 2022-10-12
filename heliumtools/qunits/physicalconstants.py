# -*- coding: utf-8 -*-

"""Constants from scipy.

See https://docs.scipy.org/doc/scipy/reference/constants.html

"""

try:
    from scipy import constants as const

    have_scipy = True
except ImportError:
    have_scipy = False


class PhysConst(object):
    """Namespace for physical constants"""

    available_constants = {
        "Avogadro constant": "na",
        "Bohr magneton": "mu_b",
        "Bohr radius": "a_b",
        "Boltzmann constant": "kb",
        "Planck constant": "h",
        "Planck constant over 2 pi": "hbar",
        "Stefan-Boltzmann constant": "sigma",
        "atomic mass constant": "u",
        "electron g factor": "ge",
        "electron mass": "me",
        "electric constant": "eps0",
        "elementary charge": "e",
        "fine-structure constant": "alpha",
        "mag. constant": "mu0",
        "speed of light in vacuum": "c",
        "standard acceleration of gravity": "g",
    }

    fallback_values = {
        "Avogadro constant": [6.022140857e+23, "mol^-1"],
        "Bohr magneton": [9.274009994e-24, "J T^-1"],
        "Bohr radius": [5.2917721067e-11, "m"],
        "Boltzmann constant": [1.38064852e-23, "J K^-1"],
        "Planck constant": [6.62607004e-34, "J s"],
        "Planck constant over 2 pi": [1.0545718e-34, "J s"],
        "Stefan-Boltzmann constant": [5.670367e-08, "W m^-2 K^-4"],
        "atomic mass constant": [1.66053904e-27, "kg"],
        "electron g factor": [-2.00231930436182, ""],
        "electron mass": [9.10938356e-31, "kg"],
        "electric constant": [8.854187817620389e-12, "F m^-1"],
        "elementary charge": [1.6021766208e-19, "C"],
        "fine-structure constant": [0.0072973525664, ""],
        "mag. constant": [1.2566370614359173e-06, "N A^-2"],
        "speed of light in vacuum": [299792458.0, "m s^-1"],
        "standard acceleration of gravity": [9.80665, "m s^-2"],
    }

    def __init__(self, unit_namespace):
        self.units = unit_namespace  # do not call it u, this is the atomic mass unit
        if have_scipy:
            for name, shortname in self.available_constants.items():
                val, unitstring, _ = const.physical_constants[name]
                self._addconstant(shortname, val, unitstring)
        else:
            for name, shortname in self.available_constants.items():
                val, unitstring = self.fallback_values[name]
                self._addconstant(shortname, val, unitstring)

    def _addconstant(self, name, val, unitstring):
        unit = self.units.from_string(unitstring)
        setattr(self, name, val * unit)


if __name__ == "__main__":
    import unitnamespace

    u = unitnamespace.UnitNamespace("si")
    const = PhysConst(u)
    # print(const.available_constants)
    for name, shortname in const.available_constants.items():
        print("{n}: {v}".format(n=name, v=getattr(const, shortname)))
