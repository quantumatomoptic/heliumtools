# -*- coding: utf-8 -*-
# cython: profile=True

"""Define the Quantity classes in cython.

Note that the interplay between numpy ndarrays, ScalarQuantity and ArrayQuantity is a
delicate thing that is controlled via the __array_priority__ attributes, NotImplemented
returns and ArrayQuantities __array_ufunc__.

"""

cimport cython
from cpython.array cimport array, copy

import numpy as np
cimport numpy as np


class EIncompatibleUnits(Exception):
    pass


class ESignatureAlreadyRegistered(Exception):
    pass


cdef class ScalarQuantity


ctypedef double[7] SiArray


cdef list symbols = ['m', 'kg', 's', 'A', 'K', 'cd', 'mol']


# ------ Helper functions to manage Quantities ------


cdef inline int is_scalarquantity(var):
    """Checks whether var is an instance of type ScalarQuantity."""
    return isinstance(var, ScalarQuantity)


@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline ScalarQuantity assert_scalarquantity(x):
    """Check if x is a ScalarQuantity, otherwise cast x as one."""
    if is_scalarquantity(x):
        return x
    else:
        return ScalarQuantity.__new__(ScalarQuantity, x)


@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline ScalarQuantity norm_quant_from_si(si_tuple):
    """Create a new ScalarQuantity with unity magnitude from a given `si_tuple`.

    The `si_tuple` holds the exponent of the si quantities representing the unit of the
    quantity in the following order: ['m', 'kg', 's', 'A', 'K', 'cd', 'mol']

    """
    ans = ScalarQuantity.__new__(ScalarQuantity, 1.0)
    ans.set_si_representation(si_tuple)
    return ans


@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline ScalarQuantity mul_scalar(x, y):
    """Helper function to multiply the ScalarQuantities x and y.

    """
    cdef ScalarQuantity xq
    cdef ScalarQuantity yq
    cdef ScalarQuantity ans
    cdef int i
    xq = assert_scalarquantity(x)
    yq = assert_scalarquantity(y)
    ans = ScalarQuantity.__new__(ScalarQuantity, xq.magnitude * yq.magnitude)
    for i from 0 <= i < 7:
        ans.si_representation[i] = xq.si_representation[i] + yq.si_representation[i]
    return ans


def mul_array(arr, x):
    """Helper function to multiply an array `arr` and a ScalarQuantity `x`.

    """
    cdef ScalarQuantity xq
    xq = assert_scalarquantity(x)
    return ArrayQuantity(arr, xq)

def div_array(arr, x):
    """Helper function to divide an array `arr` by a ScalarQuantity `x`.

    """
    cdef ScalarQuantity xq
    xq = assert_scalarquantity(x)
    return ArrayQuantity(arr, pow_scalar(xq, -1))

def div_array2(arr, x):
    """Helper function to divide a ScalarQuantity `x` by an array `arr`.

    """
    cdef ScalarQuantity xq
    xq = assert_scalarquantity(x)
    return ArrayQuantity(1/arr, xq)


@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline ScalarQuantity div_scalar(x, y):
    """Helper function to divide the ScalarQuantities x and y.

    """
    cdef ScalarQuantity xq
    cdef ScalarQuantity yq
    cdef ScalarQuantity ans
    cdef int i
    xq = assert_scalarquantity(x)
    yq = assert_scalarquantity(y)
    ans = ScalarQuantity.__new__(ScalarQuantity, xq.magnitude / yq.magnitude)
    for i from 0 <= i < 7:
        ans.si_representation[i] = xq.si_representation[i] - yq.si_representation[i]
    return ans


@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline ScalarQuantity pow_scalar(x, y):
    """Helper function to calculate x**y.

    """
    cdef ScalarQuantity xq = assert_scalarquantity(x)
    cdef ScalarQuantity yq = assert_scalarquantity(y)
    if not is_dimensionless(yq):
        raise EIncompatibleUnits('The exponent must be dimensionless.')
    cdef ScalarQuantity ans = ScalarQuantity.__new__(ScalarQuantity, xq.magnitude ** yq.magnitude)
    copy_si_representation(xq, ans, yq.magnitude)
    return ans


cdef inline void copy_si_representation(ScalarQuantity source, ScalarQuantity dest, float power):
    """Copy the si_representation of the unit of ScalarQuantity `source` to the
    ScalarQuantity `dest` while multiplying the exponents by `power`.

    """
    cdef int i
    for i from 0 <= i < 7:
        # FIXME: Find a better solution than brute force rounding. Problem is that
        # without it (u.m**3)**(1/3) != u.m
        dest.si_representation[i] = np.round(source.si_representation[i] * power, 5)


cdef inline int same_si_representation(ScalarQuantity self, ScalarQuantity other) except -1:
    """Check if two ScalarQuantities have the same si_representation, i.e. unit.

    Raise `EIncompatibleUnits` error if the units are not the same.

    """
    cdef int i
    for i from 0 <= i < 7:
        if self.si_representation[i] != other.si_representation[i]:
            raise EIncompatibleUnits(
                'Incompatible units: {} and {}'.format(self, other))
    return 0


cdef inline int is_dimensionless(ScalarQuantity x):
    """Check if the ScalarQuantity `x` is dimensionless.

    """
    cdef int i
    res = True
    for i from 0 <= i < 7:
        if x.si_representation[i] != 0:
            res = False
    return res


# ------ Global dictionaries ------

# register different unit types. All units with the same si representation belong to
# the same category.
UNITCATEGORIES = {}
cpdef add_unitcategory(si_representation, str name):
    """Add a new unit category identified by the `si_representation` and identified by `name`.

    """
    if tuple(si_representation) in UNITCATEGORIES:
        raise ESignatureAlreadyRegistered('The si representation {} is already registered as {}'
                .format(str(tuple(si_representation)), UNITCATEGORIES[tuple(si_representation)]))
    UNITCATEGORIES[tuple(si_representation)] = name


# Representation cache that links a si_representation to a representation in a compatible
# unit.
REPRESENTCACHE = {}
def setrepresent(quantity, as_unit=None, symbol='', convert_function=None, format_spec='.4g'):
    """By default, the target representation is arrived by dividing
    the current unit magnitude by the target unit magnitude, and
    appending the desired representation symbol.

    However, if a conversion_function is supplied, then INSTEAD the
    conversion function will be called so:

        output_magnitude = conversion_function(self.magnitude)
        output_symbol = symbol

        result = '{} {}'.format(output_magnitude, output_symbol)

    The intention of the function argument is to allow
    non-proportional conversion, typically temperature but also things
    like barg, bara, etc.

    Note that if a convert_function is supplied, the as_unit arg
    is IGNORED.

    """
    if not (as_unit or convert_function):
        raise Exception('Either a target unit or a conversion function must be supplied.')

    if convert_function == None:
        def proportional_conversion(instance, _):
            return instance.to(as_unit)
        convert_function = proportional_conversion
    REPRESENTCACHE[tuple(quantity.to_list())] = dict(
        convert_function=convert_function,
        symbol=symbol,
        format_spec=format_spec)


# The unit registry is a lookup list where you can find a specific Unit
# (a ScalarQuantity) from a particular symbol (a string).  Multiple entries in the
# UNITREGISTRY can point to the same unit, because there can be many synonyms for a
# particular unit, e.g. s, sec, secs, seconds
UNITREGISTRY = {}


# ------ Helper functions to manage Quantities ------

@cython.freelist(8)
cdef class ScalarQuantity:
    cdef readonly double magnitude
    cdef SiArray si_representation

    # The following is black magic to force numpy not to handle ufuncs when a ScalarQuantity
    # is involved. See https://stackoverflow.com/questions/38229953/array-and-rmul-operator-in-python-numpy
    # This in particular solves the problem that nparray * Scalarquantity results
    # in a ndarray of quantities and not a ArrayQuantity.
    __array_priority__ = 25

    def __cinit__(self, magnitude):
        self.magnitude = magnitude
        self.si_representation[:] = [0,0,0,0,0,0,0]

    def to_list(self):
        """Get the si_representation of the ScalarQuantity.

        """
        cdef list out
        cdef int i
        out = [0.0]*7
        for i from 0 <= i < 7:
            out[i] = self.si_representation[i]
        return out

    cdef inline tuple _si_representation_tuple(self):
        """Get the si_representation of the ScalarQuantity as a tuple.

        """
        return tuple(self.to_list())

    def set_si_representation(self, list si_repr):
        """Set the si_representation of the ScalarQuantity.

        """
        cdef int i
        for i from 0 <= i < 7:
            self.si_representation[i] = si_repr[i]

    def to(self, ScalarQuantity target_unit):
        """Express the ScalarQuantity in another compatible unit and strip the units.

        """
        if not is_scalarquantity(target_unit):
            raise TypeError('Target must be a unit.')
        same_si_representation(self, target_unit)
        return self.magnitude / target_unit.magnitude

    def assert_unit(self, unit):
        """Make sure the unit is compatible with the given unit.

        Otherwise rise EIncompatibleUnits error.
        """
        same_si_representation(self, unit)

    @property
    def unit(self):
        """Get the unit."""
        return self._new_quantity(1)

    def getsymbol(self):
        """Get the symbol representing the ScalarQuantity. The symbol is obtained either
        from the REPRESENTCACHE, or build from scatch.

        """
        if self._si_representation_tuple() in REPRESENTCACHE:
            r = REPRESENTCACHE[self._si_representation_tuple()]
            ret = '{}'.format(r['symbol'])
            return ret
        else:
            text = ' '.join(['{}^{}'.format(k,v) for k, v in zip(symbols, self.to_list()) if v != 0])
            ret = '{}'.format(text)
            return ret

    def _getmagnitude_represent(self):
        """Get the magnitude of the ScalarQuantity. Use the convertfunction of the
        REPRESENTCACHE if an entry for the ScalarQuantity exists.

        """
        if self._si_representation_tuple() in REPRESENTCACHE:
            r = REPRESENTCACHE[self._si_representation_tuple()]
            return r['convert_function'](self, self.magnitude)

        else:
            return self.magnitude

    def getrepresenttuple(self):
        """Get the full represent tuple of the ScalarQuantity. This includes the magnitude
        in the registered representation, the symbol and the format specification
        of for the magnitude

        """
        if self._si_representation_tuple() in REPRESENTCACHE:
            r = REPRESENTCACHE[self._si_representation_tuple()]
            format_spec = r['format_spec']
        else:
            format_spec = ''
        return self._getmagnitude_represent(), self.getsymbol(), format_spec

    def unitcategory(self):
        """Return the category of the unit."""
        if self._si_representation_tuple() in UNITCATEGORIES:
            return UNITCATEGORIES[self._si_representation_tuple()]
        else:
            msg = 'The collection of units: "{}" has not been defined as a category yet.'
            raise Exception(msg.format(str(self._si_representation_tuple())))

    def copy(self):
        """Create a copy of the ScalarQuantity.

        """
        return self._new_quantity(self.magnitude)

    def _new_quantity(self, magnitude):
        cdef ScalarQuantity ans = ScalarQuantity.__new__(ScalarQuantity, magnitude)
        copy_si_representation(self, ans, 1)
        return ans

    # ------ Define the necessary __ methods. ------

    def __str__(self):
        mag, symbol, format_spec = self.getrepresenttuple()
        number_part = format(mag, format_spec)
        if symbol == '':
            return number_part
        else:
            return ' '.join([number_part, symbol])

    def __repr__(self):
        return str(self)

    def __format__(self, format_spec):
        # Ignore the stored format_spec, use the given one.
        mag, symbol, stored_format_spec = self.getrepresenttuple()
        if format_spec == '':
            format_spec = stored_format_spec
        number_part = format(mag, format_spec)
        if symbol == '':
            return number_part
        else:
            return ' '.join([number_part, symbol])

    def __float__(self):
        if not is_dimensionless(self):
            raise EIncompatibleUnits('Must be dimensionless for __float__()')
        return self.magnitude

    def __int__(self):
        if not is_dimensionless(self):
            raise EIncompatibleUnits('Must be dimensionless for __int__()')
        return int(self.magnitude)

    def __round__(self, ndigits=None):
        return self._new_quantity(np.round(self.magnitude, ndigits))

    def __rshift__(self, other):
        """This implements the >> operator.

        Notes
        -----
        Example:
            4*u.m >> u.cm
            400

        """
        return self.to(other)

    def __reduce__(self):
        """Pickling support

        see https://docs.python.org/3/library/pickle.html#object.__reduce__
        """
        return (ScalarQuantity, (self.magnitude,), self.to_list(), None, None)

    def __setstate__(self, si_repr):
        """Pickling support

        see https://docs.python.org/3/library/pickle.html#object.__setstate__
        """
        self.set_si_representation(si_repr)

    def __hash__(self):
        return hash((self.magnitude, self._si_representation_tuple()))

    # Arithmetric for standard python types.
    # See https://cython.readthedocs.io/en/latest/src/userguide/special_methods.html

    def __add__(x, y):
        if isinstance(y, ArrayQuantity):
            # Returning NotImplemented here means python will call the __add__ function
            # of y. The same applies for the other NonImplemented returns below.
            return NotImplemented
        cdef ScalarQuantity xq
        cdef ScalarQuantity yq
        cdef ScalarQuantity ans

        xq = assert_scalarquantity(x)
        yq = assert_scalarquantity(y)
        ans = ScalarQuantity.__new__(ScalarQuantity, xq.magnitude + yq.magnitude)
        same_si_representation(xq, yq)
        copy_si_representation(xq, ans, 1)
        return ans

    def __sub__(x, y):
        if isinstance(y, ArrayQuantity):
            return NotImplemented
        cdef ScalarQuantity xq
        cdef ScalarQuantity yq
        cdef ScalarQuantity ans

        xq = assert_scalarquantity(x)
        yq = assert_scalarquantity(y)
        ans = ScalarQuantity.__new__(ScalarQuantity, xq.magnitude - yq.magnitude)
        same_si_representation(xq, yq)
        copy_si_representation(xq, ans, 1)
        return ans

    def __mul__(x, y):
        if (not isinstance(x, np.ndarray)) and (not isinstance(y, np.ndarray)):
            # scalar multiplication
            return mul_scalar(x, y)
        elif isinstance(y, ArrayQuantity):
            return NotImplemented
        elif isinstance(y, np.ndarray):
            return mul_array(y, x)
        elif isinstance(x, ArrayQuantity):
            return NotImplemented
        elif isinstance(x, np.ndarray):
            return mul_array(x, y)
        else:
            raise Exception('IMPOSSIBLE')

    def __div__(x, y):
        if (not isinstance(x, np.ndarray)) and (not isinstance(y, np.ndarray)):
            return div_scalar(x, y)
        elif isinstance(y, ArrayQuantity):
            return NotImplemented
        elif isinstance(y, np.ndarray):
            return div_array2(y, x)
        elif isinstance(x, ArrayQuantity):
            return NotImplemented
        elif isinstance(x, np.ndarray):
            return div_array(x, y)
        else:
            raise Exception('IMPOSSIBLE')

    def __truediv__(x, y):
        if (not isinstance(x, np.ndarray)) and (not isinstance(y, np.ndarray)):
            return div_scalar(x, y)
        elif isinstance(y, ArrayQuantity):
            return NotImplemented
        elif isinstance(y, np.ndarray):
            return div_array2(y, x)
        elif isinstance(x, ArrayQuantity):
            return NotImplemented
        elif isinstance(x, np.ndarray):
            return div_array(x, y)
        else:
            raise Exception('IMPOSSIBLE')

    def __pow__(x, y, z):
        return pow_scalar(x, y)

    def __neg__(x):
        cdef ScalarQuantity ans = ScalarQuantity.__new__(ScalarQuantity, -x.magnitude)
        copy_si_representation(x, ans, 1)
        return ans

    def __pos__(x):
        cdef ScalarQuantity ans = ScalarQuantity.__new__(ScalarQuantity, x.magnitude)
        copy_si_representation(x, ans, 1)
        return ans

    def __richcmp__(x, y, int op):
        """
        <   0
        <=  1
        ==  2
        !=  3
        >   4
        >=  5
        """
        if isinstance(y, ArrayQuantity):
            return NotImplemented
        cdef ScalarQuantity xq = assert_scalarquantity(x)
        cdef ScalarQuantity yq = assert_scalarquantity(y)
        same_si_representation(xq, yq)
        if op == 0:
            return xq.magnitude < yq.magnitude
        elif op == 1:
            return xq.magnitude <= yq.magnitude
        elif op == 2:
            return xq.magnitude == yq.magnitude
        elif op == 3:
            return xq.magnitude != yq.magnitude
        elif op == 4:
            return xq.magnitude > yq.magnitude
        elif op == 5:
            return xq.magnitude >= yq.magnitude

    def __abs__(x):
        cdef ScalarQuantity ans = ScalarQuantity.__new__(ScalarQuantity, abs(x.magnitude))
        copy_si_representation(x, ans, 1)
        return ans

    # ----- Numpy compatibility ------

    # TODO: For best compatibility, one should implement __array_ufunc__

    def sqrt(self):
        cdef ScalarQuantity ans = ScalarQuantity.__new__(ScalarQuantity, np.sqrt(self.magnitude))
        copy_si_representation(self, ans, 0.5)
        return ans

    def square(self):
        cdef ScalarQuantity ans = ScalarQuantity.__new__(ScalarQuantity, np.square(self.magnitude))
        copy_si_representation(self, ans, 2)
        return ans

    def cbrt(self):
        cdef ScalarQuantity ans = ScalarQuantity.__new__(ScalarQuantity, np.cbrt(self.magnitude))
        copy_si_representation(self, ans, 1./3.)
        return ans

    def reciprocal(self):
        cdef ScalarQuantity ans = ScalarQuantity.__new__(ScalarQuantity, np.reciprocal(self.magnitude))
        copy_si_representation(self, ans, -1)
        return ans

    def sin(self):
        if not is_dimensionless(self):
            raise EIncompatibleUnits('Must be dimensionless.')
        cdef ScalarQuantity ans = ScalarQuantity.__new__(ScalarQuantity, np.sin(self.magnitude))
        return ans

    def cos(self):
        if not is_dimensionless(self):
            raise EIncompatibleUnits('Must be dimensionless.')
        cdef ScalarQuantity ans = ScalarQuantity.__new__(ScalarQuantity, np.cos(self.magnitude))
        return ans

    def tan(self):
        if not is_dimensionless(self):
            raise EIncompatibleUnits('Must be dimensionless.')
        cdef ScalarQuantity ans = ScalarQuantity.__new__(ScalarQuantity, np.tan(self.magnitude))
        return ans

    def arcsin(self):
        if not is_dimensionless(self):
            raise EIncompatibleUnits('Must be dimensionless.')
        cdef ScalarQuantity ans = ScalarQuantity.__new__(ScalarQuantity, np.arcsin(self.magnitude))
        return ans

    def arccos(self):
        if not is_dimensionless(self):
            raise EIncompatibleUnits('Must be dimensionless.')
        cdef ScalarQuantity ans = ScalarQuantity.__new__(ScalarQuantity, np.arccos(self.magnitude))
        return ans

    def arctan(self):
        if not is_dimensionless(self):
            raise EIncompatibleUnits('Must be dimensionless.')
        cdef ScalarQuantity ans = ScalarQuantity.__new__(ScalarQuantity, np.arctan(self.magnitude))
        return ans

    def sinh(self):
        if not is_dimensionless(self):
            raise EIncompatibleUnits('Must be dimensionless.')
        cdef ScalarQuantity ans = ScalarQuantity.__new__(ScalarQuantity, np.sinh(self.magnitude))
        return ans

    def cosh(self):
        if not is_dimensionless(self):
            raise EIncompatibleUnits('Must be dimensionless.')
        cdef ScalarQuantity ans = ScalarQuantity.__new__(ScalarQuantity, np.cosh(self.magnitude))
        return ans

    def tanh(self):
        if not is_dimensionless(self):
            raise EIncompatibleUnits('Must be dimensionless.')
        cdef ScalarQuantity ans = ScalarQuantity.__new__(ScalarQuantity, np.tanh(self.magnitude))
        return ans

    def arcsinh(self):
        if not is_dimensionless(self):
            raise EIncompatibleUnits('Must be dimensionless.')
        cdef ScalarQuantity ans = ScalarQuantity.__new__(ScalarQuantity, np.arcsinh(self.magnitude))
        return ans

    def arccosh(self):
        if not is_dimensionless(self):
            raise EIncompatibleUnits('Must be dimensionless.')
        cdef ScalarQuantity ans = ScalarQuantity.__new__(ScalarQuantity, np.arccosh(self.magnitude))
        return ans

    def arctanh(self):
        if not is_dimensionless(self):
            raise EIncompatibleUnits('Must be dimensionless.')
        cdef ScalarQuantity ans = ScalarQuantity.__new__(ScalarQuantity, np.arctanh(self.magnitude))
        return ans

    def exp(self):
        if not is_dimensionless(self):
            raise EIncompatibleUnits('Must be dimensionless.')
        cdef ScalarQuantity ans = ScalarQuantity.__new__(ScalarQuantity, np.exp(self.magnitude))
        return ans

    def log(self):
        if not is_dimensionless(self):
            raise EIncompatibleUnits('Must be dimensionless.')
        cdef ScalarQuantity ans = ScalarQuantity.__new__(ScalarQuantity, np.log(self.magnitude))
        return ans

    def log2(self):
        if not is_dimensionless(self):
            raise EIncompatibleUnits('Must be dimensionless.')
        cdef ScalarQuantity ans = ScalarQuantity.__new__(ScalarQuantity, np.log2(self.magnitude))
        return ans

    def log10(self):
        if not is_dimensionless(self):
            raise EIncompatibleUnits('Must be dimensionless.')
        cdef ScalarQuantity ans = ScalarQuantity.__new__(ScalarQuantity, np.log10(self.magnitude))
        return ans

    def deg2rad(self):
        if not is_dimensionless(self):
            raise EIncompatibleUnits('Must be dimensionless.')
        cdef ScalarQuantity ans = ScalarQuantity.__new__(ScalarQuantity, np.deg2rad(self.magnitude))
        return ans

    def rad2deg(self):
        if not is_dimensionless(self):
            raise EIncompatibleUnits('Must be dimensionless.')
        cdef ScalarQuantity ans = ScalarQuantity.__new__(ScalarQuantity, np.rad2deg(self.magnitude))
        return ans


# dict to define the dimension of a ufunc result. The tuple specifies first the
# relation of the inputs, than of the outputs
# for all ufuncs see https://docs.scipy.org/doc/numpy-1.15.0/reference/ufuncs.html
ufunc_dims = {
        'add' : ('same', 'asfirst'),
        'subtract' : ('same', 'asfirst'),
        'multiply' : ('any', 'mult'),
        'divide' : ('any', 'div'),
        'logaddexp' : ('dimensionless','asfirst'),
        'logaddexp2' : ('dimensionless','asfirst'),
        'true_divide' : ('any', 'div'),
        'floor_divide' : ('any', 'div'),
        'negative' : ('any', 'asfirst'),
        'positive' : ('any', 'asfirst'),
        'power' : ('power', 'power'),
        'remainder' : ('any', 'asfirst'),
        'mod' : ('any', 'asfirst'),
        'fmod' : ('any', 'asfirst'),
        'divmod' : ('any', 'divmod'),
        'absolute' : ('any', 'asfirst'),
        'fabs' : ('any', 'asfirst'),
        'rint' : ('any', 'asfirst'),
        'sign' : ('any', 'none'),
        'heaviside' : ('any', 'none'),
        'conj' : ('dimensionless', 'asfirst'),
        'exp' : ('dimensionless', 'asfirst'),
        'exp2' : ('dimensionless', 'asfirst'),
        'log' : ('dimensionless', 'asfirst'),
        'log2' : ('dimensionless', 'asfirst'),
        'log10' : ('dimensionless', 'asfirst'),
        'expm1' : ('dimensionless', 'asfirst'),
        'log1p' : ('dimensionless', 'asfirst'),
        'sqrt' : ('any', 'sqrt'),
        'square' : ('any', 'square'),
        'cbrt' : ('any', 'cbrt'),
        'reciprocal' : ('any', 'reciprocal'),
        'gcd' : ('same', 'asfirst'),
        'lcm' : ('same', 'asfirst'),
        'sin' : ('dimensionless', 'asfirst'),
        'cos' : ('dimensionless', 'asfirst'),
        'tan' : ('dimensionless', 'asfirst'),
        'arcsin' : ('dimensionless', 'asfirst'),
        'arccos' : ('dimensionless', 'asfirst'),
        'arctan' : ('dimensionless', 'asfirst'),
        'arctan2' : ('same', 'dimensionless'),
        'hypot' : ('same', 'dimensionless'),
        'sinh' : ('dimensionless', 'asfirst'),
        'cosh' : ('dimensionless', 'asfirst'),
        'tanh' : ('dimensionless', 'asfirst'),
        'arcsinh' : ('dimensionless', 'asfirst'),
        'arccosh' : ('dimensionless', 'asfirst'),
        'arctanh' : ('dimensionless', 'asfirst'),
        'deg2rad' : ('dimensionless', 'asfirst'),
        'rad2deg' : ('dimensionless', 'asfirst'),
        'bitwise_and' : ('dimensionless', 'none'),
        'bitwise_or' : ('dimensionless', 'none'),
        'bitwise_xor' : ('dimensionless', 'none'),
        'invert' : ('dimensionless', 'none'),
        'left_shift' : ('dimensionless', 'none'),
        'right_shift' : ('dimensionless', 'none'),
        'greater' : ('same', 'none'),
        'greater_equal' : ('same', 'none'),
        'less' : ('same', 'none'),
        'less_equal' : ('same', 'none'),
        'not_equal' : ('same', 'none'),
        'equal' : ('same', 'none'),
        'logical_and' : ('same', 'none'),
        'logical_or' : ('same', 'none'),
        'logical_xor' : ('same', 'none'),
        'logical_not' : ('same', 'none'),
        'maximum' : ('any', 'asfirst'),
        'minimum' : ('any', 'asfirst'),
        'fmax' : ('any', 'asfirst'),
        'fmin' : ('any', 'asfirst'),
        'isfinite' : ('any', 'none'),
        'isinf' : ('any', 'none'),
        'isnan' : ('any', 'none'),
        'isnat' : ('any', 'none'),
        'fabs' : ('any', 'asfirst'),
        'signbit' : ('any', 'none'),
        'copysign' : ('any', 'asfirst'),
        'nextafter' : ('same', 'asfirst'),
        'spacing' : ('any', 'asfirst'),
        'modf' : ('any', 'asfirst'),
        'ldexp' : ('ldexp', 'asfirst'),
        'frexp' : ('any', 'frexp'),
        'floor' : ('any', 'asfirst'),
        'ceil' : ('any', 'asfirst'),
        'trunc' : ('any', 'asfirst')
        }


class ArrayQuantity(np.ndarray):
    """A quantity type for numpy arrays.

    see https://docs.scipy.org/doc/numpy-1.13.0/user/basics.subclassing.html
    """

    # Similar black magic as for ScalarQuantity. Here the array priority must be larger
    # than the one of the ScalarQuantity
    __array_priority__ = 100

    def __new__(cls, input_array, unit):
        """In the numpy machinery this is only called when when a new instance is created
        via python standard instance creation.

        cls is the ArrayQuantity class.

        """
        # magnitude of unit in si representation
        if unit is not None:
            magn = unit.magnitude
        else:
            magn = 1.0
        # create an ArrayQuantity instance.
        obj = np.asarray(input_array * magn).view(cls)
        obj._setunit(unit)
        # Finally, we must return the newly created object
        return obj

    def __array_finalize__(self, obj):
        """In the numpy machinery, this function is always called at instantiation.

        The value of obj depends on the circumstances.
        - obj is None if an instance is created following Python standard instantiation.
        - obj is an ndarray or any other ndarray subclass if the instance is created at
          a call of `array_instance.view(ArrayQuantity)`.
        - obj is an ArrayQuantity if a new instance is created in slicing.

        """
        # print('In __array_finalize__:')
        # print('   obj is %s' % repr(obj))
        # print('   type of obj is %s' % repr(type(obj)))
        # Standard instatiation. The unit is already set in new.
        if obj is None: return
        # Slicing. Set the unit.
        if isinstance(obj, ArrayQuantity):
            self._setunit(obj.unit)
        else:
            pass

    @property
    def unit(self):
        """Get the unit."""
        return getattr(self, '_unit', None)

    def _setunit(self, unit):
        """Set the unit (internal use only).

        """
        if unit is not None:
            if not is_scalarquantity(unit):
                raise TypeError('Unit must be a ScalarQuantity.')
            # si representation of unit
            si_repr = unit.to_list()
            normalized_unit = ScalarQuantity(1)
            normalized_unit.set_si_representation(si_repr)
            self._unit = normalized_unit
        else:
            unit = None

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        """ Implement ufuncs.

        see https://docs.scipy.org/doc/numpy-1.13.0/user/basics.subclassing.html

        """
        # print('In __array_ufunc__:')
        # print('   ufunc is %s' % repr(ufunc))
        # print('   method is %s' % repr(method))
        # print('   inputs are %s' % repr(inputs))
        # print('   kwargs are %s' % repr(kwargs))

        # implement the >> (__rshift__) operator
        if ufunc.__name__ == "right_shift":
            return self.to(inputs[1])

        # check the input units
        self._checkinputs(inputs, ufunc_dims[ufunc.__name__][0])
        # get the output units
        out_units = self._getoutputunits(inputs, ufunc_dims[ufunc.__name__][1])
        # print('   out_units are %s' % repr(out_units))
        # convert inputs and outputs to ndarrays.
        args = []
        for input_ in inputs:
            if isinstance(input_, ArrayQuantity):
                args.append(input_.view(np.ndarray))
            elif isinstance(input_, ScalarQuantity):
                args.append(input_.magnitude)
            else:
                args.append(input_)

        outputs = kwargs.pop('out', None)
        if outputs:
            out_args = []
            for output in outputs:
                if isinstance(output, ArrayQuantity):
                    out_args.append(output.view(np.ndarray))
                elif isinstance(output, ScalarQuantity):
                    args.append(output.magnitude)
                else:
                    out_args.append(output)
            kwargs['out'] = tuple(out_args)
        else:
            outputs = (None,) * ufunc.nout

        # call the ufunc of ndarray
        results = super(ArrayQuantity, self).__array_ufunc__(ufunc, method, *args, **kwargs)

        if results is NotImplemented:
            return NotImplemented

        if ufunc.nout == 1:
            results = (results,)
            out_units = (out_units,)

        # attach the units to the output
        ii = 0
        res = []
        for result, output in zip(results, outputs):
            if output is not None:
                res.append(ArrayQuantity(output, out_units[ii]))
            else:
                res.append(ArrayQuantity(np.asarray(result), out_units[ii]))
            ii += 1
        res = tuple(res)

        return res[0] if len(res) == 1 else res

    def __array_wrap__(self, out_arr, context=None):
        # print('In __array_wrap__:')
        # print('   self is %s' % repr(self))
        # print('   arr is %s' % repr(out_arr))
        if context is None:
            out_arr.unit = self.unit

    def __str__(self):
        if self.unit is not None:
            mag, symbol, format_spec = self.unit.getrepresenttuple()
        else:
            symbol = ''
            mag = 1.0
        arr = self.view(np.ndarray) * mag
        arr = np.asarray(arr)
        arr_repr = np.array2string(arr, formatter={'all':lambda x: format(x, format_spec)})
        return ' '.join([arr_repr, symbol]).strip()

    def __repr__(self):
        return str(self)

    def _checkinputs(self, inputs, requirement):
        """Check the dimensions of the inputs of a ufunc."""
        inp_units = []
        for inp in inputs:
            if isinstance(inp, ArrayQuantity):
                inp_units.append(inp.unit)
            elif is_scalarquantity(inp):
                inp_units.append(inp)
            elif isinstance(inp, np.ndarray):
                inp_units.append(inp[0])
            else:
                inp_units.append(inp)
        if requirement == 'same':
            for jj in range(len(inp_units) - 1):
                q1 = assert_scalarquantity(inp_units[jj])
                q2 = assert_scalarquantity(inp_units[jj + 1])
                same_si_representation(q1, q2)
        elif requirement == 'dimensionless':
            for inp in inp_units:
                inpq = assert_scalarquantity(inp)
                is_dimensionless(inpq)
        elif requirement == 'power':
            inpq = assert_scalarquantity(inp_units[1])
            is_dimensionless(inpq)
        elif requirement == 'ldexp':
            inpq = assert_scalarquantity(inp_units[1])
            is_dimensionless(inpq)

    def _getoutputunits(self, inputs, requirement):
        """Determine the output units of a ufunc.

        Note that this function must return a normalized quantity.
        """
        inp_units = []
        for inp in inputs:
            if isinstance(inp, ArrayQuantity):
                inp_units.append(inp.unit)
            elif is_scalarquantity(inp):
                inp_units.append(inp)
            elif isinstance(inp, np.ndarray):
                inp_units.append(assert_scalarquantity(inp[0]))
            else:
                inp_units.append(assert_scalarquantity(inp))
        if requirement == 'asfirst':
            si_repr = inp_units[0].to_list()
            return norm_quant_from_si(si_repr)
        elif requirement == 'mult':
            si_repr = mul_scalar(inp_units[0], inp_units[1]).to_list()
            return norm_quant_from_si(si_repr)
        elif requirement == 'div':
            si_repr = div_scalar(inp_units[0], inp_units[1]).to_list()
            return norm_quant_from_si(si_repr)
        elif requirement == 'power':
            si_repr = pow_scalar(inp_units[0], inp_units[1]).to_list()
            return norm_quant_from_si(si_repr)
        elif requirement == 'sqrt':
            si_repr = pow_scalar(inp_units[0], 1./2.).to_list()
            return norm_quant_from_si(si_repr)
        elif requirement == 'square':
            si_repr = pow_scalar(inp_units[0], 2.).to_list()
            return norm_quant_from_si(si_repr)
        elif requirement == 'cbrt':
            si_repr = pow_scalar(inp_units[0], 2.).to_list()
            return norm_quant_from_si(si_repr)
        elif requirement == 'reciprocal':
            si_repr = pow_scalar(inp_units[0], -1.).to_list()
            return norm_quant_from_si(si_repr)
        elif requirement == 'divmod':
            si_repr = inp_units[0].to_list()
            return (norm_quant_from_si(si_repr), norm_quant_from_si(si_repr))
        elif requirement == 'frexp':
            si_repr = inp_units[0].to_list()
            return (norm_quant_from_si(si_repr), ScalarQuantity(1.0))
        elif requirement == 'dimensionless':
            return ScalarQuantity(1.0)
        elif requirement == 'none':
            return None

    @property
    def magnitude(self):
        return self.view(np.ndarray)

    def to(self, target_unit):
        if not is_scalarquantity(target_unit):
            raise TypeError('Target must be a unit.')
        same_si_representation(self.unit, target_unit)
        return self.magnitude / target_unit.magnitude

    def assert_unit(self, unit):
        self.unit.assert_unit(unit)

    def __getitem__(self, val, *args, **kwargs):
        """Make slicing work"""
        res = super(ArrayQuantity, self).__getitem__(val, *args, **kwargs)
        if isinstance(res, np.ndarray):
            return ArrayQuantity(res, self.unit)
        else:
            q = norm_quant_from_si(self.unit.to_list())
            return q * res

    def __reduce__(self):
        """Support pickling"""
        return (ArrayQuantity, (self.magnitude, self.unit), None, None, None)

