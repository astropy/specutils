# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Tests for extinction curve."""

import numpy as np
from specutils.extinction import extinction

extinction_models = ['ccm89', 'od94', 'gcc09', 'f99', 'fm07']

def test_extinction_shapes():

    for model in extinction_models:

        # single value should work
        extinction(1.e4, a_v=1., model=model)

        # multiple values should return appropriate shape
        assert extinction([1.e4], a_v=1., model=model).shape == (1,) 
        assert extinction([1.e4, 2.e4], a_v=1., model=model).shape == (2,)

# TODO: resolve discrepancy here (see notes below)
def test_extinction_ccm89():

    # U, B, V, R, I, J, H, K band effective wavelengths from CCM '89 table 3
    x_inv_microns = np.array([2.78, 2.27, 1.82, 1.43, 1.11, 0.80, 0.63, 0.46])

    # A(lambda)/A(V) for R_V = 3.1 from Table 3 of CCM '89
    ratio_true = np.array([1.569, 1.337, 1.000, 0.751, 0.479, 0.282,
                           0.190, 0.114])

    wave = 1.e4 / x_inv_microns  # wavelengths in Angstroms
    a_lambda_over_a_v = extinction(wave, a_v=1., r_v=3.1, model='ccm89')

    # So far, these are close but not exact.
    # I get: [ 1.56880904  1.32257836  1. 0.75125994  0.4780346   0.28206957
    #          0.19200814  0.11572348]

    # At the sigfigs of Table 3, the differences are:
    # [ None, 0.014, None, None, 0.001, None, 0.002, 0.002 ]
    # with B band being the most significant difference.

    # a and b can be obtained with:
    # b = extinction(wave, ebv=1., r_v=0., model='ccm89')
    # a = extinction(wave, ebv=1., r_v=1., model='ccm89') - b
    #
    # b = [ 1.90899552  0.99999783  0.         -0.36499617 -0.62299483
    #      -0.36794719 -0.25046607 -0.15095612]
    # a = [ 0.95300404  0.99999842  1.          0.86900064  0.67900067
    #       0.40076222  0.27280365  0.164419  ]
    #

    # Could be due to floating point errors in original paper?
    # Should compare to IDL routines.

