import pytest
import numpy as np

import astropy.units as u
from astropy.modeling import models
from astropy.nddata import StdDevUncertainty, NDData
from astropy.tests.helper import quantity_allclose

from ..spectra import Spectrum1D, SpectralRegion
from ..analysis import (snr_threshold)


def test_snr_threshold():

    np.random.seed(42)

    # Setup 1D spectrum
    wavelengths = np.linspace(0, 10)*u.um
    flux = 100*np.abs(np.random.randn(10))*u.Jy
    uncertainty = StdDevUncertainty(np.abs(np.random.randn(10))*u.Jy)
    spectrum = Spectrum1D(spectral_axis=wavelengths, flux=flux, uncertainty=uncertainty)

    spectrum_masked = snr_threshold(spectrum, 50)

    assert all([x==y for x,y in zip(spectrum_masked.mask, [True, False, True, True, False, False, True, True, True, False])])


    # Setup 3D spectrum
    np.random.seed(42)
    wavelengths = np.arange(0, 10)*u.um
    flux = 100*np.abs(np.random.randn(3, 4, 10))*u.Jy
    uncertainty = StdDevUncertainty(np.abs(np.random.randn(3, 4, 10))*u.Jy)
    spectrum = Spectrum1D(spectral_axis=wavelengths, flux=flux, uncertainty=uncertainty)

    spectrum_masked = snr_threshold(spectrum, 50)

    masked_true = np.array([[[ True, False, False,  True, False, False,  True,  True,  True, True],
        [False,  True, False,  True,  True, False,  True,  True,  True, True],
        [ True, False, False,  True,  True, False,  True, False,  True, True],
        [ True,  True, False,  True,  True,  True, False,  True,  True, False]],
        [[ True, False, False, False,  True,  True,  True,  True,  True, True],
        [False, False,  True,  True,  True,  True,  True, False,  True, False],
        [ True, False,  True,  True,  True,  True, False,  True, False, False],
        [ True,  True, False,  True,  True,  True, False,  True,  True, True]],
        [[ True,  True,  True, False,  True,  True,  True,  True,  True, False],
        [False,  True,  True,  True,  True,  True, False,  True, False, True],
        [ True, False, False, False, False, False,  True, False, False, False],
        [ True, False,  True,  True, False, False, False,  True,  True, True]]])

    assert all([x==y for x,y in zip(spectrum_masked.mask.ravel(), masked_true.ravel())])


    # Setup 3D NDData
    np.random.seed(42)
    flux = 100*np.abs(np.random.randn(3, 4, 10))*u.Jy
    uncertainty = StdDevUncertainty(np.abs(np.random.randn(3, 4, 10))*u.Jy)
    spectrum = NDData(data=flux, uncertainty=uncertainty)

    spectrum_masked = snr_threshold(spectrum, 50)

    masked_true = np.array([[[ True, False, False,  True, False, False,  True,  True,  True, True],
        [False,  True, False,  True,  True, False,  True,  True,  True, True],
        [ True, False, False,  True,  True, False,  True, False,  True, True],
        [ True,  True, False,  True,  True,  True, False,  True,  True, False]],
        [[ True, False, False, False,  True,  True,  True,  True,  True, True],
        [False, False,  True,  True,  True,  True,  True, False,  True, False],
        [ True, False,  True,  True,  True,  True, False,  True, False, False],
        [ True,  True, False,  True,  True,  True, False,  True,  True, True]],
        [[ True,  True,  True, False,  True,  True,  True,  True,  True, False],
        [False,  True,  True,  True,  True,  True, False,  True, False, True],
        [ True, False, False, False, False, False,  True, False, False, False],
        [ True, False,  True,  True, False, False, False,  True,  True, True]]])

    assert all([x==y for x,y in zip(spectrum_masked.mask.ravel(), masked_true.ravel())])
