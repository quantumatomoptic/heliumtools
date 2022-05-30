# -*- coding: utf-8 -*-
'''
Author   : alex
Created  : 2020-10-13 17:03:22
Modified : 2020-10-16 16:12:39

Comments : some utility functions for ODT potential calculation
'''
# -- imports
import numpy as np
import scipy.constants as csts
from numpy import pi


# -- functions

# - units and scales
def unit_str(x, prec=2, unit=''):
    if x == 0:
        return '0 %s' % unit
    # define units
    _disp_prefixes = {-18: 'a', -15: 'f', -12: 'p', -9: 'n', -6: 'µ', -3: 'm',
                      0: '', 3: 'k', 6: 'M', 9: 'G', 12: 'T'}
    # get range
    n = np.log(np.abs(x)) / np.log(10)
    k = n // 3
    if np.abs(n - 3 * (k + 1)) < 1e-10:
        k += 1
    power = np.clip(3 * k, -18, 12)
    y = x / 10**(power)
    # prepare string
    if y >= 100 and prec < 4:
        fmt = '%d %s%s'
    else:
        fmt = '%.' + str(prec) + 'g %s%s'
    out = fmt % (y, _disp_prefixes[int(power)], unit)
    return out


def unit_mult(x, unit=''):
    # define units
    _disp_prefixes = {-18: 'a', -15: 'f', -12: 'p', -9: 'n', -6: 'µ', -3: 'm',
                      0: '', 3: 'k', 6: 'M', 9: 'G', 12: 'T'}
    # get range
    x = np.abs(x)
    n = np.log(x) / np.log(10)
    k = n // 3
    if np.abs(n - 3 * (k + 1)) < 1e-10:
        k += 1
    power = np.clip(3 * k, -18, 12)
    unit_string = _disp_prefixes[int(power)] + unit
    mult = 1 / 10 ** power
    return mult, unit_string


# - 2D polynomial fit
def polyfit2D(X, Y, Z, n, print_full_res=True):
    '''
    2D polynomial fit. Extension of numpy's polyfit for 2D polynomials.

    Parameters
    ----------
    X : array_like
        x-coordinate for the data
    Y : array_like
        y-coordinate for the data
    Z : array_like
        y-coordinate for the data
    n : int
        degree of the polynomial
    print_full_res : bool, optional
        set to true to display the fit result

    Returns
    --------
    p : array_like
        coefficients for the polynomial
        (see sortp() to analyze the output)
    '''
    # -- make X & Y one-dimensionnal
    x = X.flatten()
    y = Y.flatten()
    # -- prepare the "eigenbasis" for polynomial decomposition
    A = []
    for px in range(n+1):
        for py in range(n+1):
            if px+py < n+1:
                A.append(x**px*y**py)
    # -- compute the decomposition
    A = np.array(A).T
    B = Z.flatten()
    p, r, rank, s = np.linalg.lstsq(A, B, rcond=None)
    # -- display result (if asked)
    if print_full_res:
        i = 0
        for px in range(n+1):
            for py in range(n+1):
                print('x:%i | y:%i > %.5e' % (px, py, p[i]))
                i += 1
    return p


def polyval2D(X, Y, p):
    '''
    2D polynomial values. Extension of numpy's polyval for 2D polynomials.

    Parameters
    ----------
    X : array_like
        x-coordinate for the data
    Y : array_like
        y-coordinate for the data
    p : array_like
        coefficients for the polynomial (output of polyfit2D)

    Returns
    --------
    Z : array_like
        values of the polynomial at coordinates X,Y
    '''
    u = len(p)
    n = int((-1+np.sqrt(1+8*u))/2 - 1)  # strange, but true :)
    Z = X * 0
    i = 0
    for px in range(n+1):
        for py in range(n+1):
            if px+py < n+1:
                Z += p[i]*X**px*Y**py
                i += 1
    return Z


def sortp(p):
    '''
    Sorts the result of polyfit2D according to the polynmial degrees

    Parameters
    ----------
    p : array_like
        coefficients for the polynomial (output of polyfit2D)

    Returns
    -------
    psorted : dictionnary
        contains the coefficients of the polynomial sorted by X and Y degrees.
        for instance, psorted[px][py] contains the coeffcient of the
        X ** px * Y ** py term.
    '''
    u = len(p)
    n = int((-1+np.sqrt(1+8*u))/2 - 1)
    psorted = {}
    i = 0
    for px in range(n+1):
        psorted[px] = {}
        for py in range(n+1):
            if px+py < n+1:
                psorted[px][py] = p[i]
                i += 1

    return psorted


def analyze_psort(ps, unit='µK', m=4.002602 * csts.m_u):
    '''
    Analyzes the output of sortp() to retrieve potential center and frequencies
    We consider a potential with the following form :
        U(x,y) = sum_[px, py] { ps[px][py] * x**px * y**py }
    Using the values of the ps coefficients, analyze_psort() will retrieve the
    center of the potential, its value at the center, its eigenaxes and the
    corresponding frequencies.

    Parameters
    ----------
    ps : dictionnary (output of sortp())
        contains the coefficients of the polynomial sorted by X and Y degrees.
        for instance, psorted[px][py] contains the coeffcient of the
        X ** px * Y ** py term.
    unit : str, optional
        units of the fitted potential : 'J', 'µK', 'mK', 'K'
    m : float, optional
        mass of the trapped atom (in kg)

    Returns
    -------
    result : dictionnary
        contains the results of the analysis, with explicit names
    '''

    # - find eigenaxis and frequencies

    # shorthand
    c20 = ps[2][0]
    c02 = ps[0][2]
    c11 = ps[1][1]
    # find angle and frequencies
    if c02 == c20:
        # in this special case, theta = pi / 4
        # and the derivation of omega_u/v is slightly different
        theta = pi / 4
        omega_u = np.sqrt(2 * c20 + c11)
        omega_v = np.sqrt(2 * c20 - c11)
    else:
        # general case
        theta = 0.5 * np.arctan(c11 / (c20-c02))
        omega_u = np.sqrt(2 * (c20 * np.cos(theta) ** 2
                               - c02 * np.sin(theta) ** 2)
                          / (np.cos(theta) ** 4 - np.sin(theta) ** 4))
        omega_v = np.sqrt(- 2 * (c20 * np.sin(theta) ** 2
                                 - c02 * np.cos(theta) ** 2)
                          / (np.cos(theta) ** 4 - np.sin(theta) ** 4))
    # convert to SI
    if unit == 'J':
        mult = 1 / np.sqrt(m)
    elif unit == 'K':
        mult = np.sqrt(csts.k / m)
    elif unit == 'mK':
        mult = np.sqrt(csts.k / m / 1e3)
    elif unit == 'µK':
        mult = np.sqrt(csts.k / m / 1e6)
    omega_u_SI = omega_u * mult
    omega_v_SI = omega_v * mult

    # - find center
    # corresponds to solving the system :
    # c10 = omega_v**2 * sin(theta) * v0 - omega_u**2 * cos(theta) * u0
    # c01 = -omega_v**2 * cos(theta) * v0 - omega_u**2 * sin(theta) * u0

    # shorthands
    c10 = ps[1][0]
    c01 = ps[0][1]
    # prepare system
    A = np.array([[omega_v ** 2 * np.sin(theta), -omega_u**2 * np.cos(theta)],
                 [-omega_v ** 2 * np.cos(theta), -omega_u**2 * np.sin(theta)]])
    B = np.array([c10, c01])
    # solve
    r0 = np.linalg.solve(A, B)
    v0 = r0[0]
    u0 = r0[1]
    # in x,y coordinates
    x0 = u0 * np.cos(theta) - v0 * np.sin(theta)
    y0 = u0 * np.sin(theta) + v0 * np.cos(theta)

    # - find potential value at center
    U0 = ps[0][0] - 0.5 * (omega_u ** 2 * u0 ** 2 + omega_v ** 2 * v0 ** 2)

    # - prepare output
    result = {}
    result['theta'] = theta
    result['omega_u'] = omega_u_SI
    result['omega_v'] = omega_v_SI
    result['freq_u'] = omega_u_SI / 2 / pi
    result['freq_v'] = omega_v_SI / 2 / pi
    result['x0'] = x0
    result['y0'] = y0
    result['U0'] = U0
    return result
