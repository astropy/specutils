import numpy as np

import astropy.units as u
from astropy.nddata import StdDevUncertainty
from astropy.tests.helper import assert_quantity_allclose

from ..spectra.spectrum1d import Spectrum1D

from ..manipulation.replace_segment import spline_replace


def test_replace_spline():
    wave_val = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    flux_val = np.array([2, 4, 6, 8, 10, 12, 14, 16, 18, 20])
    uncert_val = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

    input_spectrum = Spectrum1D(spectral_axis=wave_val * u.AA, flux=flux_val * u.mJy,
                                uncertainty=StdDevUncertainty(uncert_val*u.mJy))

    spline_points = [3.5, 4.7, 6.8, 7.1] * u.AA

    result = spline_replace(input_spectrum, spline_points, extrapolation_treatment='zero_fill')

    assert_quantity_allclose(result.flux, flux_val*u.mJy)
    assert_quantity_allclose(result.spectral_axis, input_spectrum.spectral_axis)
