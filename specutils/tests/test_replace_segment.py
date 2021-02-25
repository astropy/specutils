import numpy as np

import astropy.units as u
from astropy.tests.helper import assert_quantity_allclose

from ..spectra.spectrum1d import Spectrum1D
from ..manipulation.replace_segment import spline_replace


def test_replace_spline():
    wave_val = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    flux_val = np.array([2, 4, 6, 8, 10, 12, 14, 16, 18, 20])

    input_spectrum = Spectrum1D(spectral_axis=wave_val * u.AA, flux=flux_val * u.mJy)

    spline_knots = [3.5, 4.7, 6.8, 7.1] * u.AA

    # with default extrapolation, recovers the input flux
    result = spline_replace(input_spectrum, spline_knots)

    assert_quantity_allclose(result.flux, flux_val*u.mJy)
    assert_quantity_allclose(result.spectral_axis, input_spectrum.spectral_axis)

    # with zero fill extrapolation, fills with zeros.
    result = spline_replace(input_spectrum, spline_knots, extrapolation_treatment='zero_fill')

    assert_quantity_allclose(result.flux[0], 0.*u.mJy)
    assert_quantity_allclose(result.flux[1], 0.*u.mJy)
    assert_quantity_allclose(result.flux[-1], 0.*u.mJy)
    assert_quantity_allclose(result.flux[-2], 0.*u.mJy)
