#!/usr/bin/env python
# coding=utf-8

# The original misu lib before my mods performed as
#  Instantiation: 0.2-0.3 ms
#  Sum: 0.2-0.3 ms
#  Numpy array instantiation: ~80 ms
#  Numpy array sum: ~6 ms

import numpy as np
import timeit
import unitnamespace

u = unitnamespace.UnitNamespace()

a=1 * u.m
b=1 * u.cm

a_arr = np.linspace(1,100,10000) * u.m
b_arr = np.linspace(1,100,10000) * u.cm

def test_instantiation():
    q = 1 * u.m

def test_sum(a, b):
    q = a + b

def test_array():
    q = np.linspace(1,100,10000) * u.m


# call timeit
print('Instantiation: {0:2.3f} ms'.format(1000*timeit.timeit("test_instantiation()", globals=globals(), number=10000)))
print('Sum: {0:2.3f} ms'.format(1000*timeit.timeit("test_sum(a, b)", globals=globals(), number=10000)))
print('Numpy array instantiation: {0:2.3f} ms'.format(1000*timeit.timeit("test_array()", globals=globals(), number=10000)))
print('Numpy array sum: {0:2.3f} ms'.format(1000*timeit.timeit("test_sum(a_arr, b_arr)", globals=globals(), number=10000)))
