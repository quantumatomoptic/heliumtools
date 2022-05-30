#!/usr/bin/env python
# -*- mode:Python; coding: utf-8 -*-

"""
@Author: victor
@Date:   14 February 2022 @ 22:58
@Last modified by:   victor
@Last modified time: 11 March 2022 @ 15:04

Comment :
"""


# -*- coding: utf-8 -*-
"""
Author   : alex
Created  : 2020-10-13 17:25:42
Modified : 2020-10-19 09:08:18

Comments : Implements the GaussianBeam class, to compute laser intensity easily
"""
# -- imports
import numpy as np
import matplotlib.pyplot as plt
from numpy import pi


# -- functions
def intensity_gauss(r, z, w0, P=1, wavelength=1550e-9):
    """
    Computes intensity for a Gaussian beam of waist w0 and power P0 at position
    (r, z). The beam is propagating along z, and the waist is located at
    r = z = 0. Lengths should be given in meters, and powers in watts.
    Intensity is returned in W/m^2
    """
    zR = pi * w0 ** 2 / wavelength
    wz = w0 * np.sqrt(1 + z ** 2 / zR ** 2)
    I0 = 2 * P / pi / w0 ** 2
    intensity = I0 * (w0 / wz) ** 2 * np.exp(-2 * (r ** 2) / wz ** 2)
    return intensity


# -- laser class
class GaussianBeam:
    """
    A single gaussian beam
    """

    def __init__(self, **kwargs):
        """
        Object initialization, sets parameters
        """
        # -- initialize default settings
        # physical parameters
        self.wavelength = 1550e-9  # meters
        self.power = 1  # Watts
        self.waist_value = 135e-6  # meters
        # geometry
        # NB: angles correspond to canonical spherical coordinates notations
        #     default = propagationg along x
        self.waist_position = (0, 0, 0)  # x, y, z
        self.waist_shift = 0  # to define a shift along the laser axis
        self.phi = 0  # angle of propagation in x/y plane
        self.theta = pi / 2  # angle of propagation, wrt z axis
        # other
        self.label = ""
        # -- initialize object
        # update attributes based on kwargs
        self.__dict__.update(kwargs)

    def intensity(self, x, y, z):
        """
        returns intensity at point (x,y,z)
        """
        # -- shorthand
        phi = self.phi
        theta = self.theta
        x0, y0, z0 = self.waist_position

        # -- change frame
        # shift center
        xc = x - x0
        yc = y - y0
        zc = z - z0

        rc = np.array([xc, yc, zc], dtype=object)

        # unit vector of laser propagation line
        v = np.array(
            [
                np.sin(theta) * np.cos(phi),  # x
                np.sin(theta) * np.sin(phi),  # y
                np.cos(theta),
            ]
        )  # z

        # axial distance (= scalar product of rc and v)
        z_dist = np.dot(v, rc)
        # radial distance (Pythagore)
        r_dist = np.sqrt(np.abs(np.linalg.norm(rc) ** 2 - z_dist ** 2))
        # add shift to axial distance
        z_dist += self.waist_shift

        # -- return intensity
        w0 = self.waist_value
        P = self.power
        wavelength = self.wavelength
        return intensity_gauss(r_dist, z_dist, w0, P, wavelength)


# -- TESTS
if __name__ == "__main__":

    # 1D intensities
    if False:
        # -- along x
        gb = GaussianBeam()
        c = "x"
        if c == "x":
            xrange = 100e-3
            yrange = 500e-6
            zrange = 500e-6
        if c == "y":
            xrange = 500e-6
            yrange = 100e-3
            zrange = 500e-6
            gb.phi = pi / 2
        if c == "z":
            xrange = 500e-6
            yrange = 500e-6
            zrange = 100e-3
            gb.theta = 0

        x = np.linspace(-xrange, xrange, 1000)
        y = np.linspace(-yrange, yrange, 1000)
        z = np.linspace(-zrange, zrange, 1000)

        plt.figure()
        plt.subplot(131)
        plt.plot(x, gb.intensity(x, 0, 0))
        plt.subplot(132)
        plt.plot(y, gb.intensity(0, y, 0))
        plt.subplot(133)
        plt.plot(z, gb.intensity(0, 0, z))
        plt.show()

    # 2D intensities
    if True:
        # init object
        gb = GaussianBeam()
        gb.waist_value = 20e-6
        gb.phi = 0 * pi / 2
        gb.theta = pi / 4
        # plot range
        xrange = 1e-3
        yrange = 1e-3
        zrange = 1e-3
        # 1D
        x = np.linspace(-xrange, xrange, 1000)
        y = np.linspace(-yrange, yrange, 1000)
        z = np.linspace(-zrange, zrange, 1000)
        # meshgrid XY
        XYx, XYy = np.meshgrid(x, y)
        # meshgrid XZ
        XZx, XZz = np.meshgrid(x, z)
        # intensities
        XYi = gb.intensity(XYx, XYy, 0)
        XZi = gb.intensity(XZx, 0, XZz)
        # plot
        fig, ax = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)
        cax = ax[0]
        cax.pcolormesh(XYx, XYy, XYi, cmap="Spectral_r")
        cax.set_xlabel("X")
        cax.set_ylabel("Y")
        cax = ax[1]
        cax.pcolormesh(XZx, XZz, XZi, cmap="Spectral_r")
        cax.set_xlabel("X")
        cax.set_ylabel("Z")
        plt.show()
