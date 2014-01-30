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
        y = x ** 1.61
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
        a = 1.752 - 0.316*x - (0.104 / ((x-4.67)**2 + 0.341))
        b = -3.090 + 1.825*x + (1.206 / ((x-4.62)**2 + 0.263))
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
    a = 1.896 - 0.372*x - 0.0108 / ((x-4.57)**2 + 0.0422)
    b = -3.503 + 2.057*x + 0.718 / ((x-4.59)**2 + 0.0530*3.1)
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

DEF F99_X0 = 4.596
DEF F99_GAMMA = 0.99
DEF F99_C3 = 3.23
DEF F99_C4 = 0.41
DEF F99_C5 = 5.9
f99_xknots = np.array([0., 1.e4/26500., 1.e4/12200., 1.e4/6000., 1.e4/5470.,
                       1.e4/4670., 1.e4/4110., 1.e4/2700., 1.e4/2600.])

def f99k(wave, double r_v):
    cdef double c1, c2, d, x, x2, y, y2, rv2
    cdef int i, n

    c2 =  -0.824 + 4.717 / r_v
    c1 =  2.030 - 3.007 * c2

    n = wave.shape[0]
    k = np.empty(n, dtype=np.float)

    # UV: analytical function
    for i in range(0, n):
        if wave[i] < 2700.:
            x = 1.e4 / wave[i]
            x2 = x * x
            d = x2 / ((x2 - F99_X0*F99_X0)**2 + x2 * F99_GAMMA*F99_GAMMA)
            k[i] = c1 + c2 * x + F99_C3 * d
            if x >= F99_C5:
                y = x - F99_C5
                y2 = y * y
                k[i] += F99_C4 * (0.5392 * y2 + 0.05644 * y2 * y)

    # Optical/IR: spline
    from scipy.interpolate import splmake, spleval
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
        x2 = f99_xknots[i] * f99_xknots[i]
        d = x2 /((x2 - F99_X0*F99_X0)**2 + x2 * F99_GAMMA*F99_GAMMA)
        kknots[i] = c1 + c2 * f99_xknots[i] + F99_C3 * d

    spline = splmake(f99_xknots, kknots, order=3)

    mask = wave >= 2700.
    xarr = 1.e4 / wave[mask]
    k[mask] = spleval(spline, xarr)

    return k
