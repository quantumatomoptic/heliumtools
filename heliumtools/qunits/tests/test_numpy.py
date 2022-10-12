#!/usr/bin/env python
# coding=utf-8

import numpy as np
from qunits.quantities import ArrayQuantity, setrepresent
import qunits
import pickle

u = qunits.UnitNamespace('all')
setrepresent(u.kg/u.s, u.kg/u.hr, 'kg/hr')

def test_repr():
    arr = np.arange(1,3)
    a = ArrayQuantity(arr, 1*u.N)
    assert repr(a) == '[1 2] N'

def test_sum_arrays():
    arr = np.arange(1,3)
    a = ArrayQuantity(arr, 1*u.N)
    b = ArrayQuantity(arr+1, 1*u.N)
    c = a + b
    assert repr(c) == '[3 5] N'

def test_sum_2():
    a = np.arange(1,3) * u.m
    b = 3 * u.m
    c =  b + a
    d =  a + b
    assert repr(c) == '[4 5] m'
    assert repr(d) == '[4 5] m'

def test_sub():
    a = np.arange(1,3) * u.m
    b = 3 * u.m
    c =  b - a
    d =  a - b
    assert repr(c) == '[2 1] m'
    assert repr(d) == '[-2 -1] m'

def test_lmult():
    arr = np.arange(1,3)
    a = 2*u.kg
    b = a * arr
    assert repr(b) == '[2 4] kg'

def test_rmult():
    arr = np.arange(1,3)
    a = 2*u.kg
    b = arr * a
    assert repr(b) == '[2 4] kg'

def test_mult_2():
    arr = np.arange(1,3) * 3*u.kg
    a = 2*u.kg
    b = arr * a
    assert repr(b) == '[6.0 12.0] kg^2.0'

def test_mult_3():
    arr = np.arange(1,3) * 3*u.kg
    a = 2*u.kg
    b = a * arr
    assert repr(b) == '[6.0 12.0] kg^2.0'

def test_div_array_1():
    arr = np.arange(1,3)
    a = 2*u.kg
    b = arr / a
    assert repr(b) == '[0.5 1.0] kg^-1.0'

def test_div_array_2():
    arr = np.arange(1,3)
    a = 2*u.kg
    b = a/arr
    assert repr(b) == '[2 1] kg'

def test_div_array_3():
    arr = np.arange(1,3) * 3*u.kg
    a = 2*u.kg
    b = arr / a
    assert repr(b) == '[1.5 3.0]'

def test_power_1():
    a = 3.0
    b = np.arange(1,4) * u.m
    c = b**a
    assert repr(c) == '[1.0 8.0 27.0] m^3.0'

def test_power_2():
    a = 3.0 * u.dimensionless
    b = np.arange(1,4) * u.m
    c = b**a
    assert repr(c) == '[1.0 8.0 27.0] m^3.0'


def test_pickle_numpy():
    var = np.array([2.5, 4]) * u.kg / u.s
    pick = pickle.dumps(var)
    pick2 = var.dumps()
    res = pickle.loads(pick)
    res2 = pickle.loads(pick2)
    assert (var[0]==res[0] and var[1]==res[1])
    assert (var[0]==res2[0] and var[1]==res2[1])

def test_numpy_multiplication():
    x1 = u.kg * np.array([1, 2, 3])
    x2 = np.array([1, 2, 3]) * u.s
    x3 = 5 * u.kg
    x4 = u.m * 5
    x5 = x1 * x2
    x6 = x3 * x4
    assert repr(x1) == '[1 2 3] kg'
    assert repr(x2) == '[1 2 3] s'
    assert repr(x3) == '5 kg'
    assert repr(x4) == '5 m'
    assert repr(x5) == '[1.0 4.0 9.0] kg^1.0 s^1.0'
    assert repr(x6) == '25.0 m^1.0 kg^1.0'

def test_numpy_sum():
    x1 = 1*u.kg + np.array([1, 2, 3])*u.kg
    x2 = np.array([1, 2, 3])*u.kg + 1*u.kg
    x3 = np.array([1, 2, 3])*u.kg + np.array([1, 2, 3])*u.kg
    assert repr(x1) == '[2 3 4] kg'
    assert repr(x2) == '[2 3 4] kg'
    assert repr(x3) == '[2 4 6] kg'

def test_numpy_division():
    x1 = np.array([1, 2, 3])/(4*u.s)
    x2 = np.array([1, 2, 3])*u.kg / (2*u.kg)
    x3 = 1*u.kg / np.array([1, 2, 3])
    x4 = 1 / np.array([1, 2, 3]) / u.m
    assert repr(x1) == '[0.25 0.5 0.75] Hz'
    assert repr(x2) == '[0.5 1.0 1.5]'
    assert repr(x3) == '[1 0.5 0.3333] kg'
    print(x4)
    assert repr(x4) == '[1.0 0.5 0.3333333333333333] m^-1.0'

def test_numpy_operations():
    x = np.array([1, 2, 3]) * u.kg
    y = x / (20*u.minutes)
    z = y**2
    w = np.sqrt(y)
    print(w)
    assert repr(y) == '[3 6 9] kg/hr'
    assert repr(z) == '[6.944444444444446e-07 2.7777777777777783e-06 6.25e-06] kg^2.0 s^-2.0'
    assert repr(w) == '[0.02886751345948129 0.040824829046386304 0.05] kg^0.5 s^-0.5'

def test_numpy_sin():
    mags = np.array([ 0.08400557, 0.19897197, 0.12407021, 0.11867142])
    x = mags * u.kN
    assert np.allclose(np.sin(mags) , np.sin(x/u.kN))
    assert np.allclose(np.sin(mags), np.sin(x.to(u.kN)))

def test_numpy_addition():
    x = np.array([1, 2, 3]) * u.kg
    y = np.array([1, 2, 3]) * u.lb
    assert repr(x+y) == '[1.454 2.907 4.361] kg'
    lbval = (x+y).to(u.lb)
    assert np.allclose(lbval,
        np.array([3.20462262,  6.40924524,  9.61386787]))


def test_numpy_subtraction():
    x = np.array([1, 2, 3]) * u.kg
    y = np.array([1, 2, 3]) * u.lb
    assert repr(x-y) == '[0.5464 1.093 1.639] kg'


def test_numpy_slice():
    x = np.array([ 0.08400557, 0.19897197, 0.12407021, 0.11867142]) * u.kg/u.hr
    a = x[:2]
    b = x[3]
    assert repr(a) == '[0.08401 0.199] kg/hr'
    assert repr(b) == '0.1187 kg/hr'

def test_abs():
    mags = np.array([ -1.2, 1.6]) * u.kN
    assert repr(abs(mags)) == '[1200 1600] N'
    assert repr(np.abs(mags)) == '[1200 1600] N'

def test_numpy_reshape():
    a = np.linspace(1,10,8) * u.mV
    b = np.reshape(a, (2,2,2))
    print(repr(b))
    assert repr(b) == '[[[0.001 0.002286]\n  [0.003571 0.004857]]\n\n [[0.006143 0.007429]\n  [0.008714 0.01]]] V'

def test_numpy_sum():
    x = np.asarray([1, 2, 3])
    y = x * u.V
    assert repr(np.sum(y)) == '6 V'

def test_numpy_var():
    a = np.linspace(1,2,8) * u.mol
    b = np.var(a)
    assert repr(b) == '0.10714285714285716 mol^2.0'

def test_numpy_std():
    a = np.linspace(1,2,8) * u.cd
    b = np.std(a)
    assert repr(b) == '0.3273 cd'

def test_zero_mult():
    """docstring for test_stuff"""
    a = np.asarray([4]) * u.m
    b = 0.0 * a
    assert repr(b) == '[0] m'

def test_sqrt_scalar():
    a = 4 * u.m**2
    b = np.sqrt(a)
    assert repr(b) == '2 m'

def test_abs_scalar():
    a = -4 * u.m
    b = np.abs(a)
    assert repr(b) == '4 m'

def test_sin_scalar():
    a = 0.5 * u.dimensionless
    b = np.sin(a)
    assert repr(b) == '0.479425538604203'

if __name__ == '__main__':
    pass
