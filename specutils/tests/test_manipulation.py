import operator
import pytest
import numpy as np

import astropy.units as u
from astropy.modeling import models
from astropy.nddata import StdDevUncertainty, NDData
from astropy.tests.helper import quantity_allclose
from specutils.wcs.wcs_wrapper import WCSWrapper

from ..spectra import Spectrum1D, SpectralRegion, SpectrumCollection
from ..manipulation import snr_threshold


def test_snr_threshold():

    np.random.seed(42)

    # Setup 1D spectrum
    wavelengths = np.arange(0, 10)*u.um
    flux = 100*np.abs(np.random.randn(10))*u.Jy
    uncertainty = StdDevUncertainty(np.abs(np.random.randn(10))*u.Jy)
    spectrum = Spectrum1D(spectral_axis=wavelengths, flux=flux, uncertainty=uncertainty)

    spectrum_masked = snr_threshold(spectrum, 50)
    assert all([x==y for x,y in zip(spectrum_masked.mask, [True, False, True, True, False, False, True, True, True, False])])

    spectrum_masked = snr_threshold(spectrum, 50, operator.gt)
    assert all([x==y for x,y in zip(spectrum_masked.mask, [True, False, True, True, False, False, True, True, True, False])])

    spectrum_masked = snr_threshold(spectrum, 50, '>')
    assert all([x==y for x,y in zip(spectrum_masked.mask, [True, False, True, True, False, False, True, True, True, False])])

    spectrum_masked = snr_threshold(spectrum, 50, operator.ge)
    assert all([x==y for x,y in zip(spectrum_masked.mask, [True, False, True, True, False, False, True, True, True, False])])

    spectrum_masked = snr_threshold(spectrum, 50, '>=')
    assert all([x==y for x,y in zip(spectrum_masked.mask, [True, False, True, True, False, False, True, True, True, False])])

    spectrum_masked = snr_threshold(spectrum, 50, operator.lt)
    assert all([not x==y for x,y in zip(spectrum_masked.mask, [True, False, True, True, False, False, True, True, True, False])])

    spectrum_masked = snr_threshold(spectrum, 50, '<')
    assert all([not x==y for x,y in zip(spectrum_masked.mask, [True, False, True, True, False, False, True, True, True, False])])

    spectrum_masked = snr_threshold(spectrum, 50, operator.le)
    assert all([not x==y for x,y in zip(spectrum_masked.mask, [True, False, True, True, False, False, True, True, True, False])])

    spectrum_masked = snr_threshold(spectrum, 50, '<=')
    assert all([not x==y for x,y in zip(spectrum_masked.mask, [True, False, True, True, False, False, True, True, True, False])])

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


    # Test SpectralCollection
    np.random.seed(42)
    flux = u.Quantity(np.random.sample((5, 10)), unit='Jy')
    spectral_axis = u.Quantity(np.arange(50).reshape((5, 10)), unit='AA')
    wcs = np.array([WCSWrapper.from_array(x).wcs for x in spectral_axis])
    uncertainty = StdDevUncertainty(np.random.sample((5, 10)), unit='Jy')
    mask = np.ones((5, 10)).astype(bool)
    meta = [{'test': 5, 'info': [1, 2, 3]} for i in range(5)]

    spec_coll = SpectrumCollection(
        flux=flux, spectral_axis=spectral_axis, wcs=wcs,
        uncertainty=uncertainty, mask=mask, meta=meta)

    spec_coll_masked = snr_threshold(spec_coll, 3)
    print(spec_coll_masked.mask)

    ma = np.array([[False, False, False, False, False, False, False,  True,  True, False],
                  [False,  True, False, False, False, False, False, False,  True, False],
                  [False, False,  True, False, False, False, False,  True, False, False],
                  [False, False, False,  True,  True, False, False, False, False, False],
                  [False, False, False, False, False, False, False, False,  True, False]])

    assert all([x==y for x,y in zip(spec_coll_masked.mask.ravel(), ma.ravel())])
