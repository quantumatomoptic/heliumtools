#!/usr/bin/env python3
# -*- mode: Python; coding: utf-8 -*-

"""
@Author: alex
@Date:   19 October 2019 @ 09:06
@Last modified by:   victor
@Last modified time: 21 January 2022 @ 16:02

Comments : Implements the coils object, to compute coils magnetic fields
"""


# -- imports
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as csts
from numpy import pi
from scipy.special import ellipe, ellipk
from numpy import linspace


# -- functions
def mag_field_coil_xy(x, y, z, radius=1.0, current=1.0):
    """
    computes the magnetic field generated by a single coil located in the (x,y)
    plane, centered on (x=0, y=0, z=0).

    source :
    [1] https://mathcurve.com/courbes2d.gb/magneticcirculaire/article%20pre.pdf

    Parameters
    ----------
    x : float / array
        cartesian x-coordinate
    y : float / array
        cartesian y-coordinate
    z : float / array
        cartesian z-coordinate
    radius : float, optional
        coil radius (m)
    current : float, optional
        coil current (A)

    """
    # -- convert x, y, z to arrays
    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)

    # -- goes to spherical coordinates
    # rho and phi
    rho = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    phi = np.arctan2(y, x)
    # theta, avoiding origin
    mask = rho == 0
    rho_corrected = np.where(mask, 1, rho)
    theta = np.arccos(np.where(mask, 1, z / rho_corrected))
    # -- compute magnetic field
    # get notations from [1]
    a = radius
    r = rho * np.sin(theta)
    k = np.sqrt(4 * a * r / ((a + r) ** 2 + z ** 2))
    # compute
    A = csts.mu_0 * current / 2 / pi  # prefactor
    E = ellipe(k ** 2)  # complete elliptic integral of second kind
    K = ellipk(k ** 2)  # complete elliptic integral of first kind
    # along z
    Bz = (
        1
        / np.sqrt((a + r) ** 2 + z ** 2)
        * (K + (a ** 2 - r ** 2 - z ** 2) / ((a - r) ** 2 + z ** 2) * E)
    )
    # along r (avoiding r=0)
    z_over_r = np.divide(z, r, where=(r != 0))
    z_over_r = np.clip(z_over_r, -1e10, 1e10)
    Br = (
        z_over_r
        / np.sqrt((a + r) ** 2 + z ** 2)
        * (-K + (a ** 2 + r ** 2 + z ** 2) / ((a - r) ** 2 + z ** 2) * E)
    )

    # -- convert to x, y
    Bx = Br * np.cos(phi)
    By = Br * np.sin(phi)

    # -- in Tesla
    Bx *= A
    By *= A
    Bz *= A

    return (Bx, By, Bz)


# -- coils objects
class SingleCoil:
    """
    A single coil, in a given plane (XY, XY or XZ)
    """

    def __init__(self, **kwargs):
        """
        Object initialization, sets parameters
        """
        # -- initialize default settings
        # physical parameters
        self.n_turns = 100  # number of turns
        self.current = 1  # coil current (A)
        # geometry
        # NB: we consider that the coil lies in a given canonical plane
        #     (XY, XZ or YZ). The coil is centered - for instance, for a coil
        #     located in the XY plane, its center is located at (x=0, y=0) -
        #     and can be shifted in the direction orthogonal to the plane
        #     - for a coil in the XY plane, the position of the coil along the
        #     Z axis can be set using the 'axial_shift' parameter
        self.radius = 10e-2  # coil radius (m)
        self.axial_shift = 30e-2  # axial shift (m)
        self.plane = "YZ"  # coil plane
        self.radial_shift = (0, 0)  # the first number match the the first letter
        #   in the plane while the second define the second letter
        # other
        self.label = ""
        self.ellipsoidal_shape = False  # if true, we must add a sqrt(2) factor when
        # computing the field (this is usefull for the IoffePritchard trap class we
        # use on the helium experiment @LCF)

        # -- initialize object
        # update attributes based on kwargs
        self.__dict__.update(kwargs)

    def field(self, x, y, z, unit="T"):
        """
        returns magnetic field at point (x,y,z)
        """

        # -- analyze inputs
        msg = "Invalid coil plane : should be 'XY', 'XZ' or 'YZ'"
        assert self.plane.upper() in ["XY", "XZ", "YZ"], msg

        unit_conversion = {"T": 1, "G": 1e4}
        unit_available = " or ".join(list(unit_conversion.keys()))
        msg = "Invalid unit choice : should be {}".format(unit_available)
        assert unit in unit_available, msg

        # -- define good coordinates
        if self.plane.upper() == "XY":
            u = x - self.radial_shift[0]
            v = y - self.radial_shift[1]
            w = z - self.axial_shift
        elif self.plane.upper() == "XZ":
            u = x - self.radial_shift[0]
            v = z - self.radial_shift[1]
            w = -y + self.axial_shift
        elif self.plane.upper() == "YZ":
            u = y - self.radial_shift[0]
            v = z - self.radial_shift[1]
            w = x - self.axial_shift

        # -- compute (on turn)
        Bu, Bv, Bw = mag_field_coil_xy(
            u, v, w, radius=self.radius, current=self.current
        )

        # -- back to good coords
        if self.plane.upper() == "XY":
            Bx = Bu
            By = Bv
            Bz = Bw
        elif self.plane.upper() == "XZ":
            Bx = Bu
            Bz = Bv
            By = -Bw
        elif self.plane.upper() == "YZ":
            By = Bu
            Bz = Bv
            Bx = Bw
        # -- good units + all turns
        Bx *= unit_conversion[unit] * self.n_turns
        By *= unit_conversion[unit] * self.n_turns
        Bz *= unit_conversion[unit] * self.n_turns

        return Bx, By, Bz


class CoilSet:
    """
    A set of coils. Designed to be handled as a SingleCoil in a transparent way
    """

    def __init__(self, coils_settings=[], label=""):
        """
        Object initialization, sets parameters

        Parameters
        ----------
        coils_settings : list, optional
            List of dictionnaries, containing settings for the collection
            of SingleCoil() constituting the CoilSet()
        label : str, optional
            Label of the CoilSet

        """
        # -- initialize settings
        # physical parameters
        self.label = label  # number of turns
        self.coil_list = []  # contains the coil collection

        # -- populate coil list
        self.update_list(coils_settings)

    def reset(self):
        self.coil_list = []

    def update_list(self, coils_settings):
        for cs in coils_settings:
            coil = SingleCoil(**cs)
            self.coil_list.append(coil)

    def field(self, x, y, z, **kwargs):
        """
        returns magnetic field at point (x,y,z)
        """
        # initialize magentic fields
        Bx = np.zeros_like(x * y * z, dtype=float)
        By = np.zeros_like(x * y * z, dtype=float)
        Bz = np.zeros_like(x * y * z, dtype=float)

        # populate
        for coil in self.coil_list:
            bx, by, bz = coil.field(x, y, z, **kwargs)
            Bx += bx
            By += by
            Bz += bz

        return Bx, By, Bz


class Coils:
    """
    Parameters
    ----------
    config : soit
    --> "MOT" seules les bobines de compensation sont traversées par un
    courant opposé (configuration anti-Helmotz)
    --> "MT" toutes les bobines sont traversées par I_top et celles du bas par
    I_bottom
    --> "Bias" seules les bobines de compensation sont traversées par un courant
    dans le même sens (configuration Helmotz).
    Axes
    ----------
    L'axe Z correspond à la gravité, et est orienté vers le haut.
    L'axe X correspond à l'axe pompe et est orienté vers le mur.
    L'axe Y correspond à l'axe Zeeman et est orienté vers les alims.
    Cela signifie que toutes les bobines sont dans le plan YZ.
    Le cluster 1 est du côté Bureau : il est donc à des x négatifs.
    Le cluster 1 est côté bureau, le 2 côté mur.
    """

    def __init__(self, **kwargs):
        ##First : default geometry for the config
        self._config = "MOT"  # possible choice : MOT, MT, Bias
        self._i_top = 200
        self._i_bottom = 0
        #  thickness of the copper wires
        self.wire_thickness = 4e-3
        # distance between the coils and the center of the vaccuum chamber
        self.cluster_distance = 26e-3
        # distance from the X axis of the quadrupoles.
        self.quadrupole_axe_distance = 45e-3
        # default current, in ampere
        self.currents = {
            "compensation": 200,
            "dipole": 0,
            "quadrupole": 0,
        }
        self.diameters = {
            "compensation": 68 * 2e-3,  # maximum diameter
            "dipole": 23 * 2e-3,  # minimum diameter
            "quadrupole": 40e-3,  # this is the maximal diameter of the elliptical quadrupoles
        }
        self._center_distance = {
            "compensation": 54e-3,  # expérimental ?
            "dipole": 50e-3,  # minimum diameter
            "quadrupole": 26e-3,
        }
        self._coils = {
            "compensation": 0,
            "dipole": 0,
            "quadrupole": 0,
        }  # this is a dictionnary of CoilSet (for now, it is 0)
        # update the default parameters
        self._coils_properties = {
            "compensation": 0,
            "dipole": 0,
            "quadrupole": 0,
        }  # here we store the coils properties before we turn them into a Coil Set
        self.__dict__.update(kwargs)

        # Now generate the coils cluster

    def update_currents(self):
        if self.config in ["MOT", "Bias"]:
            self.currents["compensation"] = self._i_top
            self.currents["dipole"] = 0
            self.currents["quadrupole"] = 0
        elif self.config == "MT":
            self.currents["compensation"] = self._i_top
            self.currents["dipole"] = self._i_top + self._i_bottom
            self.currents["quadrupole"] = self._i_top + self._i_bottom

    def update(self):
        self.update_compensation_coils()
        self.update_dipole_coils()
        self.update_quadrupole_coils()

    @property
    def compensation(self):
        self.update_compensation_coils()
        return self._coils["compensation"]

    @property
    def dipole(self):
        self.update_dipole_coils()
        return self._coils["dipole"]

    @property
    def quadrupole(self):
        self.update_quadrupole_coils()
        return self._coils["quadrupole"]

    @property
    def cluster(self):
        self.update_compensation_coils()
        self.update_dipole_coils()
        self.update_quadrupole_coils()
        quadru = self._coils_properties["quadrupole"]
        comp = self._coils_properties["compensation"]
        dip = self._coils_properties["dipole"]
        return CoilSet(quadru + comp + dip, label="Cluster")

    @property
    def coils(self):
        self.update()
        coils_list = []
        for name, coil in self._coils_properties.items():
            coils_list.append(coil)
        return CoilSet(coils_list, label="Every coils")

    @property
    def config(self):
        return self._config

    @property
    def i_top(self):
        return self._i_top

    @property
    def i_bottom(self):
        return self._i_bottom

    def set_config(self, value):
        # check the config required
        msg = "Invalid configuration. The configuration you require ({}) does not exist. It should be 'MOT', 'MT' or 'Bias'.".format(
            value
        )
        assert value in ["MOT", "MT", "Bias"], msg
        self._config = value

    def set_i_top(self, value):
        self._i_top = value

    def set_i_bottom(self, value):
        self._i_bottom = value

    def update_compensation_coils(self):
        """
        updates the compensation coils cluster. It must be called each time a current is changed.
        """
        self.update_currents()
        Ic = -self.currents["compensation"]  # --> le courant est opposé au dipole
        e = self.wire_thickness  # copper wire thickness
        L = self._center_distance["compensation"]
        d_comp = self.diameters["compensation"]  # inner diameter of compensation coil
        Ncepes = 2
        Nclong = 6
        plane = "YZ"
        ############### START TO DEFINE COMPENSATION CLUSTER  ###############
        ## CLUSTER 1
        coil = [
            {
                "plane": plane,
                "n_turns": 1,
                "current": Ic,
                "radius": (d_comp - (2 * it2) * e) / 2,
                "axial_shift": -L - it * e,
            }
            for it in range(Nclong)
            for it2 in range(Ncepes)
        ]
        # On change le courant dans l'autre cluster selon que la configuration
        # soit en Helmotz ou anti Helmotz.
        if self.config == "MOT":
            Ic = -Ic
        ## CLUSTER 2
        coil += [
            {
                "plane": plane,
                "n_turns": 1,
                "current": Ic,
                "radius": (d_comp - (2 * it2) * e) / 2,
                "axial_shift": L + it * e,
            }
            for it in range(Nclong)
            for it2 in range(Ncepes)
        ]
        self._coils["compensation"] = CoilSet(coil, label="Compensation")
        self._coils_properties["compensation"] = coil
        ############### END TO DEFINE COMPENSATION CLUSTER  ###############

    def update_dipole_coils(self):
        """
        updates the dipole coils cluster. It must be called each time a current is changed.
        """
        self.update_currents()
        Id = self.currents["dipole"]
        e = self.wire_thickness  # copper wire thickness
        L = self._center_distance["dipole"]
        L = 50e-3
        # quadrupoles (hence 5 thickness of copper wire)
        d_dip = self.diameters["dipole"]
        plane = "YZ"
        Ndlong = 5
        Ndepes = 4
        ## CLUSTER 1
        coil = [
            {
                "plane": plane,
                "n_turns": 1,
                "current": Id,
                "radius": (d_dip + 2 * it2 * e) / 2,
                "axial_shift": -L - it * e,
            }
            for it in range(Ndlong)
            for it2 in range(Ndepes)
        ]

        ## CLUSTER 2
        coil += [
            {
                "plane": plane,
                "n_turns": 1,
                "current": Id,
                "radius": (d_dip + (2 * it2) * e) / 2,
                "axial_shift": L + it * e,
            }
            for it in range(Ndlong)
            for it2 in range(Ndepes)
        ]
        self._coils_properties["dipole"] = coil
        self._coils["dipole"] = CoilSet(coil, label="Dipole")

    def update_quadrupole_coils(self):
        """
        updates the quadrupole coils cluster. It must be called each time a current is changed.
        Here, remimber that the coils are not circular but ellipticals so that it increases
        the value of the gradient. According to the PhD thesis of A. Browaeys, one can read
        page 108 that
        "Les valeurs du quadrupôle sont elliptiques afin d'augmenter la valeur du gradient (en
        théorie, par rapport à des bobines circulaires, on gagne un facteur sqrt(2) sur le gradient
        total en utilisant des bobines elliptiques, équivalentes à deux bobines circulaires de même
        rayon décalées de 45°)".
        To take that into account, we add a factor sqrt(2) on the current since the B fields and
        its gradient should be linear with it.
        """
        self.update_currents()
        Iq = self.currents["quadrupole"] * 1.4  # this is the root(2) factor.
        e = self.wire_thickness  # copper wire thickness
        L = self._center_distance["quadrupole"]
        r0_quad = self.quadrupole_axe_distance  # distance of the quadrupole center
        # to the X axis
        d_quad = self.diameters["quadrupole"]
        plane = "YZ"
        Nqlong = 5
        Nqepes = 3
        ############### START TO DEFINE QUADRUPOLE CLUSTER  ###############
        # fmt: off

        coilQ1A = [{"plane": plane, "n_turns": 1,"current": Iq, "radius": (d_quad - ( 2 * it2) * e) / 2, "axial_shift": L +  it * e, "radial_shift": (r0_quad, 0)} for it in range(Nqlong) for it2 in range(Nqepes) ]
        coilQ2A = [{"plane": plane, "n_turns": 1,"current": -Iq, "radius": (d_quad - (2 * it2) * e) / 2, "axial_shift": -L  - it * e, "radial_shift": (r0_quad, 0)} for it in range(Nqlong) for it2 in range(Nqepes) ]
        coilQ1C = [{"plane": plane, "n_turns": 1,"current": Iq, "radius": (d_quad - (2 * it2) * e) / 2, "axial_shift": L + + it * e, "radial_shift": (-r0_quad, 0)} for it in range(Nqlong) for it2 in range(Nqepes) ]
        coilQ2C = [{"plane": plane, "n_turns": 1,"current": -Iq, "radius": (d_quad - ( 2 * it2) * e) / 2, "axial_shift": -L - it * e, "radial_shift": (-r0_quad, 0)} for it in range(Nqlong) for it2 in range(Nqepes) ]
        coilQ1B = [{"plane": plane, "n_turns": 1,"current": -Iq, "radius": (d_quad - (2 * it2) * e) / 2, "axial_shift": L  + it * e, "radial_shift": (0,r0_quad)} for it in range(Nqlong) for it2 in range(Nqepes) ]
        coilQ2B = [{"plane": plane, "n_turns": 1,"current": Iq, "radius": (d_quad - ( 2 * it2) * e) / 2, "axial_shift": -L - it * e, "radial_shift": (0,r0_quad)} for it in range(Nqlong) for it2 in range(Nqepes) ]
        coilQ1D = [{"plane": plane, "n_turns": 1,"current": -Iq, "radius": (d_quad -(2 * it2) * e) / 2, "axial_shift": L + it * e, "radial_shift": (0,-r0_quad)} for it in range(Nqlong) for it2 in range(Nqepes) ]
        coilQ2D = [{"plane": plane, "n_turns": 1,"current": Iq, "radius": (d_quad - ( 2 * it2) * e) / 2, "axial_shift": -L  - it * e, "radial_shift": (0,-r0_quad)} for it in range(Nqlong) for it2 in range(Nqepes) ]
        # coilQ1A = coilQ1A[0:-1] --> pour prendre en compte l'endomagement d'une bobine ??
        # fmt: on
        ############### END TO DEFINE QUADRUPOLE CLUSTER  ###############
        self._coils["quadrupole"] = CoilSet(
            coilQ1A
            + coilQ1B
            + coilQ1C
            + coilQ1D
            + coilQ2A
            + coilQ2B
            + coilQ2C
            + coilQ2D,
            label="Quadrupole",
        )
        self._coils_properties["quadrupole"] = (
            coilQ1A
            + coilQ1B
            + coilQ1C
            + coilQ1D
            + coilQ2A
            + coilQ2B
            + coilQ2C
            + coilQ2D
        )


# -- TESTS
if __name__ == "__main__":
    coils = Coils()
    if False:
        # Here you can see the system
        magpy.displaySystem(coils.compensation)

    ##########
    if False:
        coils.set_config("MT")  # on règle la configuration des IGBT
        coils.set_i_top(210)
        coils.set_i_bottom(40)
        cluster = coils.cluster

        ## Affichage
        NX = 100
        NY = 101
        NZ = 100
        xmax = 15
        xs = np.linspace(-xmax, xmax, NX)
        ys = np.linspace(-xmax, xmax, NY)
        zs = np.linspace(-xmax, xmax, NZ)

        ##Plots 1D
        fig, ax = plt.subplots(3, figsize=(15, 8))
        Bx, By, Bz = cluster.field(0.001 * xs, 0, 0, unit="G")
        B = np.sqrt(Bx ** 2 + Bz ** 2 + By ** 2)
        print(np.min(B))
        ax[0].plot(xs, B)
        ax[0].set_title("Selon X")
        # gradient
        gradB = 0.5 * 10 * (B[2:] - B[0:-2]) / (xs[1] - xs[0])
        ax2x = ax[0].twinx()
        ax2x.set_ylabel("Grad(B) (G/cm)")
        ax2x.plot(
            xs[1:-1],
            gradB,
            "*g",
        )
        ## Y ##
        Bx, By, Bz = cluster.field(0, 0.001 * ys, 0, unit="G")
        B = np.sqrt(Bx ** 2 + Bz ** 2 + By ** 2)
        print(np.min(B))
        ax[1].plot(ys, B)
        ax[1].set_title("Selon Y")
        # gradient
        gradB = 0.5 * 10 * (B[2:] - B[0:-2]) / (ys[1] - ys[0])
        ax2y = ax[1].twinx()
        ax2y.set_ylabel("Grad(B) (G/cm)")
        ax2y.plot(
            ys[1:-1],
            gradB,
            "*g",
        )
        ## Z ##
        Bx, By, Bz = cluster.field(0, 0, 0.001 * zs, unit="G")
        B = np.sqrt(Bx ** 2 + Bz ** 2 + By ** 2)
        ax[2].plot(zs, B)
        ax[2].set_title("Selon Z")
        # gradient
        print(np.min(B))
        gradB = 0.5 * 10 * (B[2:] - B[0:-2]) / (zs[1] - zs[0])
        ax2z = ax[2].twinx()
        ax2z.set_ylabel("Grad(B) (G/cm)")
        ax2z.plot(
            zs[1:-1],
            gradB,
            "*g",
        )

        for a in ax.flat:
            a.set(xlabel="Position (mm)", ylabel="|B| (G)")
        for a in ax.flat:
            a.label_outer()
        plt.show()
    if True:
        coils.set_config("MT")  # on règle la configuration des IGBT
        coils.set_i_top(200)
        coils.set_i_bottom(50)
        cluster = coils.dipole

        ## Affichage
        NX = 100
        NY = 100
        NZ = 100
        xmax = 40
        xs = np.linspace(-xmax, xmax, NX)
        ys = np.linspace(-xmax, xmax, NY)
        zs = np.linspace(-xmax, xmax, NZ)
        ##Plots 2D
        fig, axs = plt.subplots(1, 2, figsize=(15, 6))
        x, y = np.meshgrid(xs, ys)
        Bx, By, Bz = cluster.field(0.001 * x, 0.001 * y, 0, unit="G")
        B = np.sqrt(Bx ** 2 + Bz ** 2 + By ** 2)
        axs[0].streamplot(x, y, Bx, By, color="w")
        pcm = axs[0].pcolormesh(x, y, B, vmin=0, cmap="Spectral_r")
        axs[0].set_ylabel("Y (mm)")
        axs[0].set_xlabel("X (mm)")
        plt.colorbar(pcm, ax=axs[0], label="|B| (Gauss)")

        y, z = np.meshgrid(ys, zs)
        Bx, By, Bz = cluster.field(0, 0.001 * y, 0.001 * z, unit="G")
        B = np.sqrt(Bx ** 2 + Bz ** 2 + By ** 2)
        y, z = np.meshgrid(ys, zs)
        axs[1].streamplot(y, z, By, Bz, color="w")
        pcm = axs[1].pcolormesh(y, z, B, vmin=0, cmap="Spectral_r")
        axs[1].set_ylabel("Z (mm)")
        axs[1].set_xlabel("Y (mm)")
        plt.colorbar(pcm, ax=axs[1], label="|B| (Gauss)")
        plt.show()

    if True:
        comp_coils = coils.compensation
        print(comp_coils)
        # - grids and computation
        # grid
        x = np.linspace(-2, 2, 500)
        z = np.linspace(-2, 2, 500)
        y = 0
        x, z = np.meshgrid(x, z)
        # compute
        Bx, By, Bz = comp_coils.field(0.001 * x, y, 0.001 * z, unit="G")
        B = np.sqrt(Bx ** 2 + Bz ** 2 + By ** 2)

        # print
        print("|B|_max = %.2f G" % np.max(B))
        for [ax, Bi] in zip(["x", "y", "z"], [Bx, By, Bz]):
            print("|B%s|_max = %.2f G" % (ax, np.max(Bi)))

        # plot
        fig, ax = plt.subplots(1, 1, figsize=(10, 8))
        ax.streamplot(x, z, Bx, Bz, color="k")
        pcm = ax.pcolormesh(
            x, z, np.sqrt(Bx ** 2 + Bz ** 2 + By ** 2), vmin=0, cmap="Spectral_r"
        )
        ax.set_ylabel("Z (m)")
        ax.set_xlabel("X (m)")
        plt.colorbar(pcm, ax=ax, label="|B| (Gauss)")
        plt.show()
    if False:
        coils.set_config("MT")  # on règle la configuration des IGBT
        coils.set_i_top(198)
        coils.set_i_bottom(-20)
        cluster = coils.cluster

        ## Affichage
        NX = 100
        NY = 101
        NZ = 100
        xmax = 1
        xs = np.linspace(-xmax, xmax, NX)
        ys = np.linspace(-xmax, xmax, NY)
        zs = np.linspace(-xmax, xmax, NZ)

        ##Plots 1D selon Y
        fig, ax = plt.subplots(4, figsize=(15, 8))
        Bx, By, Bz = cluster.field(0.001, 0.001 * ys, 0, unit="G")
        B = np.sqrt(Bx ** 2 + Bz ** 2 + By ** 2)
        print(np.min(B))
        ax2cluster = ax[0].twinx()
        ax2cluster.set_ylabel("B (G)")
        ax2cluster.plot(ys, Bx, "seagreen", label="Bx", linestyle="-")
        ax2cluster.plot(ys, By, "orange", label="By", linestyle="--")
        ax2cluster.plot(ys, Bz, "darkred", label="Bz", linestyle="-.")
        ax2cluster.legend()
        ax[0].plot(ys, B)
        ax[0].set_title("Cluster")

        Bx, By, Bz = coils.dipole.field(0, 0.001 * ys, 0, unit="G")
        B = np.sqrt(Bx ** 2 + Bz ** 2 + By ** 2)
        ax[1].plot(ys, B)
        ax[1].set_title("Dipoles")
        ax2dipole = ax[1].twinx()
        ax2dipole.set_ylabel("B (G)")
        ax2dipole.plot(ys, Bx, "seagreen", label="Bx", linestyle="-")
        ax2dipole.plot(ys, By, "orange", label="By", linestyle="--")
        ax2dipole.plot(ys, Bz, "darkred", label="Bz", linestyle="-.")

        Bx, By, Bz = coils.quadrupole.field(0, 0.001 * ys, 0, unit="G")
        B = np.sqrt(Bx ** 2 + Bz ** 2 + By ** 2)
        ax[2].plot(ys, B)
        ax[2].set_title("Quadrupoles")
        ax2quadru = ax[2].twinx()
        ax2quadru.set_ylabel("B (G)")
        ax2quadru.plot(ys, Bx, "seagreen", label="Bx", linestyle="-")
        ax2quadru.plot(ys, By, "orange", label="By", linestyle="--")
        ax2quadru.plot(ys, Bz, "darkred", label="Bz", linestyle="-.")

        Bx, By, Bz = coils.compensation.field(0, 0.001 * ys, 0, unit="G")
        B = np.sqrt(Bx ** 2 + Bz ** 2 + By ** 2)
        ax[3].plot(ys, B)
        ax[3].set_title("Compensation")
        ax2comp = ax[3].twinx()
        ax2comp.set_ylabel("B (G)")
        ax2comp.plot(ys, Bx, "seagreen", label="Bx", linestyle="-")
        ax2comp.plot(ys, By, "orange", label="By", linestyle="--")
        ax2comp.plot(ys, Bz, "darkred", label="Bz", linestyle="-.")

        for a in ax.flat:
            a.set(xlabel="Position sur l'axe Y (mm)", ylabel="|B| (G)")
        for a in ax.flat:
            a.label_outer()
        plt.show()
