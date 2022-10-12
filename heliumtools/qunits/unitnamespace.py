#!/usr/bin/env python
# coding=utf-8

"""Define the UnitNameSpace.

"""

from __future__ import division, print_function
from collections import OrderedDict  # OrderedDict: from py 2.7 on
import sys
import re
import functools
import json
import inspect
import os.path as osp

from heliumtools.qunits.siprefixes import siprefixes_sym
import heliumtools.qunits.quantities
from heliumtools.qunits.quantities import ScalarQuantity, ESignatureAlreadyRegistered


class UnitNamespace(object):
    """A namespace for all defined units."""

    si_symbols = ["m", "kg", "s", "A", "K", "cd", "mol"]

    def __init__(self, context="si"):
        """Initialize the unit namespace.

        We create all SI base units plus the dimensionless unit.
        See https://physics.nist.gov/cuu/Units/units.html

        Parameters
        ----------
        context : string (default: 'si')
            In which context are we working? Used to restrict the available units.
            The value 'all' loads all units.

        """
        # dimensionless is special
        self.dimensionless = ScalarQuantity(1.0)
        heliumtools.qunits.quantities.add_unitcategory(
            self.dimensionless.to_list(), "Dimensionless"
        )
        self.known_units = ["dimensionless"]

        # meter
        self.add_unit(
            symbols=["m", "metre", "metres", "meter", "meters"],
            si_repr=[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            scale_factor=1.0,
            unit_category="Length",
            representative_symbol="m",
            create_metric_prefixes_for=["m"],
            metric_skip_function=None,
        )
        # kg (special since prefix already included)
        self.add_unit(
            symbols=["kg", "kilograms", "kilogram"],
            si_repr=[0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            scale_factor=1.0,
            unit_category="Mass",
            representative_symbol="kg",
            create_metric_prefixes_for=[],
            metric_skip_function=None,
        )
        self.add_unit(
            symbols=["g", "grams", "gram"],
            si_repr=[0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            scale_factor=1e-3,
            create_metric_prefixes_for=["g"],
            metric_skip_function=lambda p: p in ["k"],  # kg are already there.
        )

        # seconds
        self.add_unit(
            symbols=["s", "second", "sec", "seconds", "secs"],
            si_repr=[0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
            scale_factor=1.0,
            unit_category="Time",
            representative_symbol="s",
            create_metric_prefixes_for=["s"],
            metric_skip_function=lambda p: p in ["a"],  # no "as" since it is a keyword
        )

        # ampere
        self.add_unit(
            symbols=["A", "ampere", "amperes", "amp", "amps"],
            si_repr=[0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
            scale_factor=1.0,
            unit_category="Electric current",
            representative_symbol="A",
            create_metric_prefixes_for=["A"],
            metric_skip_function=None,
        )

        # ampere
        self.add_unit(
            symbols=["K", "kelvin"],
            si_repr=[0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            scale_factor=1.0,
            unit_category="Temperature",
            representative_symbol="K",
            create_metric_prefixes_for=["K"],
            metric_skip_function=None,
        )

        # candela
        self.add_unit(
            symbols=["cd", "candela", "ca"],
            si_repr=[0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
            scale_factor=1.0,
            unit_category="Luminous intensity",
            representative_symbol="cd",
            create_metric_prefixes_for=["cd"],
            metric_skip_function=None,
        )

        # mol
        self.add_unit(
            symbols=["mol", "mole", "moles"],
            si_repr=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
            scale_factor=1.0,
            unit_category="Ammount of substance",
            representative_symbol="mol",
            create_metric_prefixes_for=["mol"],
            metric_skip_function=None,
        )

        # create derived units according to the context
        self._create_derived_units(context)

    def add_unit(
        self,
        symbols,
        si_repr,
        scale_factor,
        unit_category="",
        representative_symbol="",
        create_metric_prefixes_for=[],
        metric_skip_function=None,
    ):
        """Add a unit to the namespace.

        Parameters
        ----------
        symbols : list of unit symbols
            These will be put into the class namespace, and will be entered as keys in
            the global UNITREGISTRY.

        si_repr : list
            List containing the representation of the unit in exponents of the si units
            in the follorwing order ['m', 'kg', 's', 'A', 'K', 'cd', 'mol']

        scale_factor : float
            Scale factor **wrt the SI units** of the unit.

        unit_category : string (default: "")
            Category the unit belongs to. A unit category is defined by having a unique
            decomposition in SI base units. This should be specified for the first unit
            per category only. Otherwise a ESignatureAlreadyRegistered is triggered,
            which is caugth and translated into a warning.

        representative_symbol : string (default: "")
            Symbol that should be used to represent the unit in a result of a calculation.
            Must be in symbols. Note that the each value here overides the previously set
            representative symbol for the unit category (see unit_category above).

        create_metric_prefixes_for : list of unit symbols (default: empty list)
            List of symbols (must also be in symbols) to create derived units with the
            metric prefixes.

        metric_skip_function : callable (default: None)
            Callable that returns true for the metric prefixes for which the prefixed
            unit should not be created.

        """
        # First define the unit based on si
        unit = ScalarQuantity(scale_factor)
        unit.set_si_representation(si_repr)
        # Add to registry and namespace
        for symbol in symbols:
            symbol = symbol.strip()
            if symbol == "":
                continue
            self._add_attr_unit(symbol, unit)
            self._add_to_registry(symbol, unit)
        # Add category
        if not unit_category == "":
            # print("Adding type {} to category {}".format(unit.units(), unit_category))
            try:
                heliumtools.qunits.quantities.add_unitcategory(
                    unit.to_list(), str(unit_category)
                )
            except ESignatureAlreadyRegistered as e:
                print(
                    "WARNING: Can not resister {} for unit category {}: {}".format(
                        symbols[0], unit_category, e
                    )
                )
        # Set representative symbol
        if not representative_symbol == "":
            self._check_represent(representative_symbol, symbols)
            heliumtools.qunits.quantities.setrepresent(
                unit, as_unit=unit, symbol=representative_symbol
            )
        # Metric prefixes
        if not create_metric_prefixes_for == []:
            self._check_metric_prefix_request(create_metric_prefixes_for, symbols)
            for symbol in create_metric_prefixes_for:
                self.create_metric_prefixes(symbol, unit, metric_skip_function)

    def create_metric_prefixes(self, symbol, unit, skipfunction=None):
        """Populates the UNITREGISTRY and the namespace with all the
        SI-prefixed versions of the given symbol.

        """
        for prefix in siprefixes_sym:
            if skipfunction and skipfunction(prefix):
                continue
            prefsymb = "{p}{s}".format(p=prefix, s=symbol)
            prefunit = 10 ** (float(siprefixes_sym[prefix].exponent)) * unit
            self._add_attr_unit(prefsymb, prefunit)
            self._add_to_registry(prefsymb, prefunit)

    def _check_metric_prefix_request(self, metricpref_list, symbols):
        """Check if all the symbols we shall create metric prefixes for are in symbols."""
        if not isinstance(metricpref_list, list):
            raise TypeError(
                "The symbols for which we shall create prefixed units must be given as a list."
            )
        for mp in metricpref_list:
            if mp not in symbols:
                raise ValueError(
                    "Can't create metric prefixed units for {s}. Not in symbols {sym}".format(
                        s=mp, sym=symbols
                    )
                )

    def _check_represent(self, symb, symbols):
        """Check if the representing symbol symb is valid."""
        if symb not in symbols:
            raise ValueError(
                "Representative symbol {s} not in symbols {sym}".format(
                    s=symb, sym=symbols
                )
            )

    def _add_attr_unit(self, symb, unit):
        """Add the unit as an attribute to the namespace."""
        setattr(self, symb, unit)
        self.known_units.append(symb)

    def _add_to_registry(self, symbol, unit):
        """Add symbol representing unit to the UNITREGISTRY."""
        if symbol not in heliumtools.qunits.quantities.UNITREGISTRY.keys():
            heliumtools.qunits.quantities.UNITREGISTRY[symbol] = unit
        else:
            raise ValueError(
                "Unit symbol {s} already present in the registry.".format(s=symbol)
            )

    def _create_derived_units(self, context):
        """Create derived units according to the requested context."""
        datapath = osp.join(osp.dirname(__file__), "unitdefs.json")
        with open(datapath) as f:
            unitdefs = json.load(f, object_pairs_hook=OrderedDict)
        for unitdef in unitdefs.values():
            conts = unitdef["used in contexts"]
            if ("all" in conts) or (context in conts) or (context == "all"):
                unit = self.from_string(
                    unitdef["representation in SI or earlier defined unit"]
                )
                mag_si = self._get_si_mag(unit)
                # optional stuff
                if "category" in unitdef.keys():
                    cat = unitdef["category"]
                else:
                    cat = ""
                if "representative symbol" in unitdef.keys():
                    rep_sym = unitdef["representative symbol"]
                else:
                    rep_sym = ""
                if "metric prefixes for" in unitdef.keys():
                    prefixes = unitdef["metric prefixes for"]
                else:
                    prefixes = []
                if "skipped prefixes" in unitdef.keys():
                    if unitdef["skipped prefixes"]:
                        skipfcn = eval(
                            "lambda p: p in {}".format(unitdef["skipped prefixes"])
                        )
                    else:
                        skipfcn = None
                # add the unit
                self.add_unit(
                    symbols=unitdef["symbols"],
                    si_repr=unit.to_list(),
                    scale_factor=unitdef["scale factor"] * mag_si,
                    unit_category=cat,
                    representative_symbol=rep_sym,
                    create_metric_prefixes_for=prefixes,
                    metric_skip_function=skipfcn,
                )

    def _get_si_mag(self, unit):
        """Get the magnitude of unit when expressed in si units."""
        si_repr = unit.to_list()
        si_unit = ScalarQuantity(1)
        si_unit.set_si_representation(si_repr)
        res = unit.to(si_unit)
        return res

    def from_string(self, string):
        """Create a Unit instance from the supplied string.

        The string has to be in the format that heliumtools.qunits uses for string representations, i.e.
        the following works:

        1.0 m
        1 m
        1 m^2 s^-1
        1 m/s
        1.248e+05 m/s
        -1.158e+05 m/s kg

        """
        # empty string?
        if not string.strip() == "":
            # Multiplication: replace all whitespace surounded by a-z,A-Z,0-9 with *
            string = re.sub(r"([a-zA-Z0-9])(\s+)([a-zA-Z0-9])", r"\1*\3", string)

            # Exponentiation: replace all ^ with **
            string = re.sub(r"\^", r"**", string)

            # inject self
            string = re.sub(r"\b([a-zA-Z])", r"self.\1", string)
            # print(string)

            res = None
            try:
                res = eval(string)
            except NameError:
                raise ValueError(f"Failed to parse {string} as unit")
            except SyntaxError:
                raise ValueError(f"Failed to parse {string} as unit")
        else:
            res = self.dimensionless
        return res

    def from_list(self, si_repr):
        """Create a unit from a si representation list.

        Parameters
        ----------
        si_repr : list
            list with the units exponents in terms of si units. The order of the list is
            ['m', 'kg', 's', 'A', 'K', 'cd', 'mol']

        """
        if not isinstance(si_repr, list):
            raise TypeError("si_repr must be a list.")
        if not len(si_repr) == 7:
            raise ValueError("si_repr must have length 7.")
        unit = ScalarQuantity(1.0)
        unit.set_si_representation(si_repr)
        return unit


# ----- helpers -----


def units_to_this_ns(unit_namespace):
    """Add the units defined in unit_namespace to the callers scope.

    Use call this to have direct access to the units, i.e. in your module do:

        import heliumtools.qunits
        u = heliumtools.qunits.UnitNamespace()
        heliumtools.qunits.units_to_this_ns(u)

    Note that this uses a bit of black magic...

    """
    stack = inspect.stack()
    try:
        locals_ = stack[1][0].f_locals
    finally:
        del stack
    for symb in unit_namespace.known_units:
        locals_[symb] = getattr(unit_namespace, symb)


# ----- decorators -----


def noquantity(func):
    """Decorator to assure the input parameters are no Unit.

    Notes
    -----
    Usage example:

        @noquantity
        def example(a, b):
            # do something

    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        # call the wrapped function
        for arg in args:
            if isinstance(arg, ScalarQuantity):
                raise TypeError("Quantity arguments not allowed.")
        for arg in kwargs.values():
            if isinstance(arg, ScalarQuantity):
                raise TypeError("Quantity arguments not allowed.")
        res = func(*args, **kwargs)
        return res

    return wrapper


def calc_unitless(out_unit_list, **in_unit_kwargs):
    """Decorator to convert the input parameters to magnitude and the output
    parameter back to unit.

    Parameters
    ----------
    out_unit_list : list of Quantities
        Quantities of the output in correct order.

    in_unit_kwargs : Quantities as keyword arguments
        Use this to specify the Quantities to which the input arguments should be
        converted to.

    Notes
    -----
    Usage example:

        @dimensions([u.m/u.s, u.kg], a=u.m, b=u.s, c=u.kg)
        def example(a, b, c=u.kg):
            # do something
            return a/b, c

    """

    def calc_unitless_decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            # names of the arguments
            if sys.version_info.major == 2:
                arg_names = func.func_code.co_varnames
            elif sys.version_info.major == 3:
                arg_names = func.__code__.co_varnames
            else:
                raise Exception("Invalid Python version!")
            # transfer args to kwargs
            kwargs.update(zip(arg_names, args))

            # now do the conversion
            if not in_unit_kwargs.keys() == kwargs.keys():
                raise ValueError(
                    "Quantity keyword arguments must match the function parameter names."
                )
            conv_kwargs = {}
            for kk in in_unit_kwargs.keys():
                conv_kwargs[kk] = kwargs[kk].to(in_unit_kwargs[kk])

            # call the function and convert the result.
            res = func(**conv_kwargs)
            if isinstance(res, (list, tuple)):
                if not len(res) == len(out_unit_list):
                    raise ValueError(
                        "Number of output units must match the number of the function return values."
                    )
                return [r * out_unit_list[jj] for jj, r in enumerate(res)]
            else:
                return res * out_unit_list[0]

        return wrapper

    return calc_unitless_decorator


# This is a decorator that will ensure arguments match declared unit category
def dimensions(**_params_):
    """Decorator to assure the parameters given have the correct unit category.

    Notes
    -----
    Usage example:

        @dimensions(a='Length', b='Time')
        def example(a, b):
            # do something

    """

    def check_types(_func_, _params_=_params_):
        def modified(*args, **kw):
            if sys.version_info.major == 2:
                arg_names = _func_.func_code.co_varnames
            elif sys.version_info.major == 3:
                arg_names = _func_.__code__.co_varnames
            else:
                raise Exception("Invalid Python version!")
            kw.update(zip(arg_names, args))
            for name, category in _params_.items():
                param = kw[name]
                assert isinstance(
                    param, ScalarQuantity
                ), """Parameter "{}" must be an instance of class Quantity
(and must be of unit type "{}").""".format(
                    name, category
                )
                assert (
                    param.unitcategory() == category
                ), 'Parameter "{}" must be unit type "{}".'.format(name, category)
            return _func_(**kw)

        modified.__name__ = _func_.__name__
        modified.__doc__ = _func_.__doc__
        # Py 3 only
        # modified.__annotations__ = _func_.__annotations__
        return modified

    # For IDEs, make sure the arg lists propagate through to the user
    return check_types


# ----- helpers for units with offset -----

# We do not support units with offset. Here are some helpers for temperature values.


@noquantity
def k_val_from_c(celsius):
    kelvin = celsius - 273.15
    return kelvin


@noquantity
def c_val_from_k(kelvin):
    celsius = kelvin + 273.15
    return celsius


@noquantity
def k_val_from_f(fahrenheit):
    kelvin = (fahrenheit + 459.67) * 5 / 9
    return kelvin


@noquantity
def f_val_from_k(kelvin):
    fahrenheit = kelvin * 9 / 5 - 459.67
    return fahrenheit


@noquantity
def c_val_from_f(fahrenheit):
    celsius = (fahrenheit - 32) * 5 / 9
    return celsius


@noquantity
def f_val_from_c(celsius):
    fahrenheit = celsius * 9 / 5 + 32
    return fahrenheit


if __name__ == "__main__":
    u = UnitNamespace("all")
    units_to_this_ns(u)
    print(m)
    u.from_string("1 m^-1 s")
    print(5 * mHz)
    print(5 * kg)
    a = 5 * kg
    print(a + 3 * g >> mg)

    tt = 4 * R
    print(tt)
    heliumtools.qunits.quantities.setrepresent(K, R, "R")
    print(4 * R)
    heliumtools.qunits.quantities.setrepresent(K, K, "K")
    print(4 * R)

    k_val_from_c(5)

    @dimensions(a="Length")
    def test(a):
        return a

    test(5 * m)

    @calc_unitless([u.m / u.s, u.m], a=u.m, b=u.s)
    def test2(a, b):
        return a / b, a

    print(test2(5 * u.m, 2 * u.s))

    # information on units
    aa = u.ng
    print(aa.unitcategory())
    print(aa.unitstring())

    # test a more complex setrepresent
    bb = 1e-18 * u.J
    print(bb)
    # report energies in Hz
    heliumtools.qunits.quantities.setrepresent(
        J, symbol="Hz", convert_function=lambda u, mag: mag / 6.62607004e-34
    )
    print(bb)
