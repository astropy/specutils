import operator
import numpy as np

import astropy.units as u
from astropy.nddata import StdDevUncertainty, NDData
from ..utils.wcs_utils import gwcs_from_array

from ..spectra import Spectrum1D, SpectralRegion, SpectrumCollection
from ..manipulation import snr_threshold, excise_regions, linear_exciser


def test_true_exciser():
    np.random.seed(84)
    spectral_axis = np.linspace(5000, 5099, num=100)*u.AA
    flux = (np.random.randn(100) + 3) * u.Jy
    spec = Spectrum1D(flux=flux, spectral_axis=spectral_axis)
    region = SpectralRegion([(5005,5010), (5060,5065)]*u.AA)
    excised_spec = excise_regions(spec, region)

    assert len(excised_spec.spectral_axis) == len(spec.spectral_axis)-10
    assert len(excised_spec.flux) == len(spec.flux)-10
    assert np.isclose(excised_spec.flux.sum(), 243.2617*u.Jy, atol=0.001*u.Jy)

    excised_spec = excise_regions(spec, region, inclusive_upper=True)
    assert len(excised_spec.spectral_axis) == len(spec.spectral_axis)-12
    assert len(excised_spec.flux) == len(spec.flux)-12


def test_linear_exciser():
    np.random.seed(84)
    spectral_axis = np.linspace(5000,5099,num=100)*u.AA
    flux = (np.random.rand(100)*100) * u.Jy
    spec = Spectrum1D(flux=flux, spectral_axis = spectral_axis)
    region = SpectralRegion([(5020,5030)]*u.AA)
    excised_spec = excise_regions(spec, region, exciser = linear_exciser)

    assert len(excised_spec.spectral_axis) == len(spec.spectral_axis)
    assert len(excised_spec.flux) == len(spec.flux)
    assert np.isclose(excised_spec.flux[25], 34.9864*u.Jy, atol=0.001*u.Jy)


def test_snr_threshold():

    np.random.seed(42)

    # Setup 1D spectrum
    wavelengths = np.arange(0, 10)*u.um
    flux = 100*np.abs(np.random.randn(10))*u.Jy
    uncertainty = StdDevUncertainty(np.abs(np.random.randn(10))*u.Jy)
    spectrum = Spectrum1D(spectral_axis=wavelengths, flux=flux, uncertainty=uncertainty)

    spectrum_masked = snr_threshold(spectrum, 50)
    assert all([x==y for x,y in zip(spectrum_masked.mask, [False, True, False, False, True, True, False, False, False, True])])

    spectrum_masked = snr_threshold(spectrum, 50, operator.gt)
    assert all([x==y for x,y in zip(spectrum_masked.mask, [False, True, False, False, True, True, False, False, False, True])])

    spectrum_masked = snr_threshold(spectrum, 50, '>')
    assert all([x==y for x,y in zip(spectrum_masked.mask, [False, True, False, False, True, True, False, False, False, True])])

    spectrum_masked = snr_threshold(spectrum, 50, operator.ge)
    assert all([x==y for x,y in zip(spectrum_masked.mask, [False, True, False, False, True, True, False, False, False, True])])

    spectrum_masked = snr_threshold(spectrum, 50, '>=')
    assert all([x==y for x,y in zip(spectrum_masked.mask, [False, True, False, False, True, True, False, False, False, True])])

    spectrum_masked = snr_threshold(spectrum, 50, operator.lt)
    assert all([not x==y for x,y in zip(spectrum_masked.mask, [False, True, False, False, True, True, False, False, False, True])])

    spectrum_masked = snr_threshold(spectrum, 50, '<')
    assert all([not x==y for x,y in zip(spectrum_masked.mask, [False, True, False, False, True, True, False, False, False, True])])

    spectrum_masked = snr_threshold(spectrum, 50, operator.le)
    assert all([not x==y for x,y in zip(spectrum_masked.mask, [False, True, False, False, True, True, False, False, False, True])])

    spectrum_masked = snr_threshold(spectrum, 50, '<=')
    assert all([not x==y for x,y in zip(spectrum_masked.mask, [False, True, False, False, True, True, False, False, False, True])])

    # Setup 3D spectrum
    np.random.seed(42)
    wavelengths = np.arange(0, 10)*u.um
    flux = 100*np.abs(np.random.randn(3, 4, 10))*u.Jy
    uncertainty = StdDevUncertainty(np.abs(np.random.randn(3, 4, 10))*u.Jy)
    spectrum = Spectrum1D(spectral_axis=wavelengths, flux=flux, uncertainty=uncertainty)

    spectrum_masked = snr_threshold(spectrum, 50)

    masked_true = np.array([
        [[False, True, True, False, True, True, False, False, False, False],
         [True, False, True, False, False, True, False, False, False, False],
         [False, True, True, False, False, True, False, True, False, False],
         [False, False, True, False, False, False, True, False, False, True]],
        [[False, True, True, True, False, False, False, False, False, False],
         [True, True, False, False, False, False, False, True, False, True],
         [False, True, False, False, False, False, True, False, True, True],
         [False, False, True, False, False, False, True, False, False, False]],
        [[False, False, False, True, False, False, False, False, False, True],
         [True, False, False, False, False, False, True, False, True, False],
         [False, True, True, True, True, True, False, True, True, True],
         [False, True, False, False, True, True, True, False, False, False]]])

    assert all([x==y for x,y in zip(spectrum_masked.mask.ravel(), masked_true.ravel())])

    # Setup 3D NDData
    np.random.seed(42)
    flux = 100*np.abs(np.random.randn(3, 4, 10))*u.Jy
    uncertainty = StdDevUncertainty(np.abs(np.random.randn(3, 4, 10))*u.Jy)
    spectrum = NDData(data=flux, uncertainty=uncertainty)

    spectrum_masked = snr_threshold(spectrum, 50)

    masked_true = np.array([
        [[False, True, True, False, True, True, False, False, False, False],
         [True, False, True, False, False, True, False, False, False, False],
         [False, True, True, False, False, True, False, True, False, False],
         [False, False, True, False, False, False, True, False, False, True]],
        [[False, True, True, True, False, False, False, False, False, False],
         [True, True, False, False, False, False, False, True, False, True],
         [False, True, False, False, False, False, True, False, True, True],
         [False, False, True, False, False, False, True, False, False, False]],
        [[False, False, False, True, False, False, False, False, False, True],
         [True, False, False, False, False, False, True, False, True, False],
         [False, True, True, True, True, True, False, True, True, True],
         [False, True, False, False, True, True, True, False, False, False]]])

    assert all([x==y for x,y in zip(spectrum_masked.mask.ravel(), masked_true.ravel())])

    # Test SpectralCollection
    np.random.seed(42)
    flux = u.Quantity(np.random.sample((5, 10)), unit='Jy')
    spectral_axis = u.Quantity(np.arange(50).reshape((5, 10)), unit='AA')
    wcs = np.array([gwcs_from_array(x) for x in spectral_axis])
    uncertainty = StdDevUncertainty(np.random.sample((5, 10)), unit='Jy')
    mask = np.ones((5, 10)).astype(bool)
    meta = [{'test': 5, 'info': [1, 2, 3]} for i in range(5)]

    spec_coll = SpectrumCollection(
        flux=flux, spectral_axis=spectral_axis, wcs=wcs,
        uncertainty=uncertainty, mask=mask, meta=meta)

    spec_coll_masked = snr_threshold(spec_coll, 3)
    print(spec_coll_masked.mask)

    ma = np.array([[True, True, True, True, True, True, True, False, False, True],
                   [True, False, True, True, True, True, True, True, False, True],
                   [True, True, False, True, True, True, True, False, True, True],
                   [True, True, True, False, False, True, True, True, True, True],
                   [True, True, True, True, True, True, True, True, False, True]])

    assert all([x==y for x,y in zip(spec_coll_masked.mask.ravel(), ma.ravel())])
