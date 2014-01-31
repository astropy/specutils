# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Extinction law functions."""
import numpy as np

# Optical coefficients for CCM89-like laws:
cdef int ccm89_coeffs_n = 8
cdef double *ccm89_coeffs_a = [1., 0.17699, -0.50447, -0.02427, 0.72085,
                               0.01979, -0.77530, 0.32999]
cdef double *ccm89_coeffs_b = [0., 1.41338, 2.28305, 1.07233, -5.38434,
                               -0.62251, 5.30260, -2.09002]
cdef int od94_coeffs_n = 9
cdef double *od94_coeffs_a = [1., 0.104, -0.609, 0.701, 1.137, -1.718,
                              -0.827, 1.647, -0.505]
cdef double *od94_coeffs_b = [0., 1.952, 2.908, -3.989, -7.985, 11.102,
                              5.491, -10.805, 3.347]

cdef double ccm89like(double wave, double r_v, double optical_coeffs_a[],
                      double optical_coeffs_b[], int n):
    cdef double x, a, b, y, y2, y3, yn
    cdef int i
    x = 1.e4 / wave

    if x < 1.1:
        y = x**1.61
        a = 0.574 * y
        b = -0.527 * y
    elif x < 3.3:
        y = x - 1.82
        a = optical_coeffs_a[0]
        b = optical_coeffs_b[0]
        yn = 1.
        for i in range(1, n):
            yn *= y
            a += optical_coeffs_a[i] * yn
            b += optical_coeffs_b[i] * yn
    elif x < 8.:
        y = x - 4.67
        a = 1.752 - 0.316*x - (0.104 / (y*y + 0.341))
        y = x - 4.62
        b = -3.090 + 1.825*x + (1.206 / (y*y + 0.263))
        if x > 5.9:
            y = x - 5.9
            y2 = y * y
            y3 = y2 * y
            a += -0.04473*y2 - 0.009779*y3
            b += 0.2130*y2 + 0.1207*y3
    else:
        y = x - 8.
        y2 = y * y
        y3 = y2 * y
        a = -0.070*y3 + 0.137*y2 - 0.628*y - 1.073
        b = 0.374*y3 - 0.420*y2 + 4.257*y + 13.670

    return a + b / r_v

# UV portion of GCC09 law, used for wave < 3030.30303 (x > 3.3)
cdef double gcc09uv(double wave, double r_v):
    cdef double x, y, y2, y3, a, b
    x = 1.e4/wave
    y = x - 4.57
    a = 1.896 - 0.372*x - 0.0108 / (y*y + 0.0422)
    y = x - 4.59
    b = -3.503 + 2.057*x + 0.718 / (y*y + 0.0530*3.1)
    if x > 5.9:
        y = x - 5.9
        y2 = y * y
        y3 = y * y2
        a += -0.110 * y2 - 0.0099 * y3
        b += 0.537 * y2 + 0.0530 * y3
    return a + b / r_v

def ccm89(double[:] wave, double a_v, double r_v):
    cdef int n = wave.shape[0]
    cdef int i
    res = np.empty(n, dtype=np.float)
    for i in range(n):
        res[i] = a_v * ccm89like(wave[i], r_v, ccm89_coeffs_a,
                                 ccm89_coeffs_b, ccm89_coeffs_n)
    return res

def od94(double[:] wave, double a_v, double r_v):
    cdef int n = wave.shape[0]
    cdef int i
    res = np.empty(n, dtype=np.float)
    for i in range(n):
        res[i] = a_v * ccm89like(wave[i], r_v, od94_coeffs_a, od94_coeffs_b,
                                 od94_coeffs_n)
    return res

def gcc09(double[:] wave, double a_v, double r_v):
    cdef int n = wave.shape[0]
    cdef int i
    res = np.empty(n, dtype=np.float)
    for i in range(n):
        if wave[i] < 3030.3030303030303:
            res[i] = a_v * gcc09uv(wave[i], r_v)
        else:
            res[i] = a_v * ccm89like(wave[i], r_v, od94_coeffs_a,
                                     od94_coeffs_b, od94_coeffs_n)
    return res

# These should be the same as in extinction.py
DEF F99_X0 = 4.596
DEF F99_GAMMA = 0.99
DEF F99_C3 = 3.23
DEF F99_C4 = 0.41
DEF F99_C5 = 5.9
DEF F99_X02 = F99_X0 * F99_X0
DEF F99_GAMMA2 = F99_GAMMA * F99_GAMMA

# Used for wave < 2700.
def f99uv(double[:] wave, double a_v, double r_v):
    cdef double c1, c2, d, x, x2, y, y2, rv2, k
    cdef int i, n

    c2 =  -0.824 + 4.717 / r_v
    c1 =  2.030 - 3.007 * c2

    n = wave.shape[0]
    res = np.empty(n, dtype=np.float)
    for i in range(0, n):
        x = 1.e4 / wave[i]
        x2 = x * x
        y = x2 - F99_X02
        d = x2 / (y * y + x2 * F99_GAMMA2)
        k = c1 + c2 * x + F99_C3 * d
        if x >= F99_C5:
            y = x - F99_C5
            y2 = y * y
            k += F99_C4 * (0.5392 * y2 + 0.05644 * y2 * y)
        res[i] = a_v * (1. + k / r_v)

    return res

def f99kknots(double[:] xknots, double r_v):
    cdef double c1, c2, d, x, x2, y, rv2
    cdef int i
    c2 =  -0.824 + 4.717 / r_v
    c1 =  2.030 - 3.007 * c2
    rv2 = r_v * r_v

    kknots = np.empty(9, dtype=np.float)
    kknots[0] = -r_v
    kknots[1] = 0.26469 * r_v/3.1 - r_v
    kknots[2] = 0.82925 * r_v/3.1 - r_v
    kknots[3] = -0.422809 + 1.00270*r_v + 2.13572e-04*rv2 - r_v
    kknots[4] = -5.13540e-02 + 1.00216 * r_v - 7.35778e-05*rv2 - r_v
    kknots[5] = 0.700127 + 1.00184*r_v - 3.32598e-05*rv2 - r_v
    kknots[6] = (1.19456 + 1.01707*r_v - 5.46959e-03*rv2 +
                 7.97809e-04 * rv2 * r_v - 4.45636e-05 * rv2*rv2 - r_v)
    for i in range(7,9):
        x2 = xknots[i] * xknots[i]
        y = (x2 - F99_X02)
        d = x2 /(y * y + x2 * F99_GAMMA2)
        kknots[i] = c1 + c2*xknots[i] + F99_C3 * d

    return kknots

# These constants should be the same as in extinction.py
DEF FM07_X0 = 4.592
DEF FM07_GAMMA = 0.922
DEF FM07_C1 = -0.175
DEF FM07_C2 = 0.807
DEF FM07_C3 = 2.991
DEF FM07_C4 = 0.319
DEF FM07_C5 = 6.097
DEF FM07_X02 = FM07_X0 * FM07_X0
DEF FM07_GAMMA2 = FM07_GAMMA * FM07_GAMMA
DEF FM07_R_V = 3.1  # Fixed for the time being (used in fm07kknots)

# Used for wave < 2700.
def fm07uv(double[:] wave, double a_v):
    cdef double d, x, x2, y, k
    cdef int i, n

    n = wave.shape[0]
    res = np.empty(n, dtype=np.float)
    for i in range(0, n):
        x = 1.e4 / wave[i]
        x2 = x * x
        y = x2 - FM07_X02
        d = x2 / (y*y + x2 * FM07_GAMMA2)
        k = FM07_C1 + FM07_C2 * x + FM07_C3 * d
        if x > FM07_C5:
            y = x - FM07_C5
            k += FM07_C4 * y * y
        res[i] = a_v * (1. + k / 3.1)

    return res

# This is mainly defined here rather than as a constant in the public module
# so that we don't have to define the FM07 constants in two places.
def fm07kknots(double[:] xknots):
    cdef double d
    cdef int i, n

    n = xknots.shape[0]
    kknots = np.empty(n, dtype=np.float)
    for i in range(0, 5):
        kknots[i] = (-0.83 + 0.63*FM07_R_V) * xknots[i]**1.84 - FM07_R_V
    kknots[5] = 0.
    kknots[6] = 1.322
    kknots[7] = 2.055
    for i in range(8, 10):
        d = xknots[i]**2 / ((xknots[i]**2 - FM07_X02)**2 +
                            xknots[i]**2 * FM07_GAMMA2)
        kknots[i] = FM07_C1 + FM07_C2 * xknots[i] + FM07_C3 * d
    return kknots
