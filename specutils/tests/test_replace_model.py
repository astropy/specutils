import numpy as np

import astropy.units as u
from astropy.nddata import StdDevUncertainty
from astropy.tests.helper import assert_quantity_allclose

from ..spectra.spectrum1d import Spectrum1D
from ..manipulation.replace_model import model_replace


def test_replace_spline():
    wave_val = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    flux_val = np.array([2, 4, 6, 8, 10, 12, 14, 16, 18, 20])

    input_spectrum = Spectrum1D(spectral_axis=wave_val * u.AA, flux=flux_val * u.mJy)

    spline_knots = [3.5, 4.7, 6.8, 7.1] * u.AA

    # with default extrapolation, recovers the input flux
    result = model_replace(input_spectrum, spline_knots)

    assert result.uncertainty is None

    assert_quantity_allclose(result.flux, flux_val*u.mJy)
    assert_quantity_allclose(result.spectral_axis, input_spectrum.spectral_axis)

    # with zero fill extrapolation, fills with zeros.
    result = model_replace(input_spectrum, spline_knots, extrapolation_treatment='zero_fill')

    assert_quantity_allclose(result.flux[0], 0.*u.mJy)
    assert_quantity_allclose(result.flux[1], 0.*u.mJy)
    assert_quantity_allclose(result.flux[-1], 0.*u.mJy)
    assert_quantity_allclose(result.flux[-2], 0.*u.mJy)


def test_replace_spline_uncert():
    wave_val = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    flux_val = np.array([2, 4, 6, 8, 10, 12, 14, 16, 18, 20])
    uncert_val = flux_val / 10.

    uncert = StdDevUncertainty(uncert_val * u.mJy)

    input_spectrum = Spectrum1D(spectral_axis=wave_val * u.AA, flux=flux_val * u.mJy,
                                uncertainty=uncert)

    spline_knots = [3.5, 4.7, 6.8, 7.1] * u.AA

    result = model_replace(input_spectrum, spline_knots)

    assert isinstance(result.uncertainty, StdDevUncertainty)
    assert result.flux.unit == result.uncertainty.unit

    assert_quantity_allclose(result.uncertainty.quantity, uncert_val*u.mJy)

    # Now try with the non-default no-uncertainty mode: result should
    # have no uncertainty even when input has.
    result = model_replace(input_spectrum, spline_knots,
                           interpolate_uncertainty=False)

    assert result.uncertainty is None


def test_replace_spline_uncert_zerofill():
    wave_val = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    flux_val = np.array([2, 4, 6, 8, 10, 12, 14, 16, 18, 20])
    uncert_val = flux_val / 10.

    uncert = StdDevUncertainty(uncert_val * u.mJy)

    input_spectrum = Spectrum1D(spectral_axis=wave_val * u.AA, flux=flux_val * u.mJy,
                                uncertainty=uncert)

    spline_knots = [3.5, 4.7, 6.8, 7.1] * u.AA

    result = model_replace(input_spectrum, spline_knots,
                           extrapolation_treatment='zero_fill')

    assert isinstance(result.uncertainty, StdDevUncertainty)
    assert result.flux.unit == result.uncertainty.unit

    assert_quantity_allclose(result.uncertainty.quantity, uncert_val*u.mJy)
