import astropy.units as u
import numpy as np
import pytest

from astropy.nddata import StdDevUncertainty
from ..spectra.spectrum1d import Spectrum1D
from ..spectra.spectrum_collection import SpectrumCollection
from ..wcs.wcs_wrapper import WCSWrapper
from ..wcs.wcs_adapter import WCSAdapter


@pytest.fixture
def spectrum_collection():
    flux = u.Quantity(np.random.sample((5, 10)), unit='Jy')
    spectral_axis = u.Quantity(np.arange(50).reshape((5, 10)), unit='AA')
    wcs = np.array([WCSWrapper.from_array(x).wcs for x in spectral_axis])
    uncertainty = StdDevUncertainty(np.random.sample((5, 10)), unit='Jy')
    mask = np.ones((5, 10)).astype(bool)
    meta = [{'test': 5, 'info': [1, 2, 3]} for i in range(5)]

    spec_coll = SpectrumCollection(
        flux=flux, spectral_axis=spectral_axis, wcs=wcs,
        uncertainty=uncertainty, mask=mask, meta=meta)

    return spec_coll


def test_create_spectrum_collection(spectrum_collection):
    assert spectrum_collection.ndim == 1
    assert spectrum_collection.shape == (5, )
    assert spectrum_collection.nspectral == 10
    assert isinstance(spectrum_collection.flux, u.Quantity)
    assert isinstance(spectrum_collection.spectral_axis, u.Quantity)
    assert spectrum_collection.wavelength.unit.physical_type == 'length'
    assert spectrum_collection.frequency.unit.physical_type == 'frequency'
    assert spectrum_collection.energy.unit.physical_type == 'energy'
    assert isinstance(spectrum_collection.uncertainty, StdDevUncertainty)
    assert spectrum_collection.uncertainty.array.shape == spectrum_collection.flux.shape


def test_spectrum_collection_slicing(spectrum_collection):
    # Test that slicing a spectrum collection returns a Spectrum1D
    single_spec = spectrum_collection[0]

    assert isinstance(single_spec, Spectrum1D)
    assert isinstance(single_spec.flux, u.Quantity)
    assert isinstance(single_spec.spectral_axis, u.Quantity)
    assert single_spec.flux.shape == spectrum_collection.flux.shape[1:]
    assert isinstance(single_spec.wcs, WCSAdapter)
    assert single_spec.mask.shape == single_spec.flux.shape
    assert len(single_spec.meta) == len(spectrum_collection.meta[0])
    assert isinstance(single_spec.uncertainty, StdDevUncertainty)


def test_collection_without_optional_arguments():
    # Without uncertainties
    flux = u.Quantity(np.random.sample((5, 10)), unit='Jy')
    spectral_axis = u.Quantity(np.arange(50).reshape((5, 10)), unit='AA')
    uncertainty = StdDevUncertainty(np.random.sample((5, 10)), unit='Jy')
    wcs = np.array([WCSWrapper.from_array(x).wcs for x in spectral_axis])
    mask = np.ones((5, 10)).astype(bool)
    meta = [{'test': 5, 'info': [1, 2, 3]} for i in range(5)]

    spec_coll = SpectrumCollection(
        flux=flux, spectral_axis=spectral_axis, wcs=wcs, mask=mask, meta=meta)

    # Without mask
    spec_coll = SpectrumCollection(
        flux=flux, spectral_axis=spectral_axis, uncertainty=uncertainty, 
        wcs=wcs, meta=meta)

    # Without meta
    spec_coll = SpectrumCollection(
        flux=flux, spectral_axis=spectral_axis, uncertainty=uncertainty, 
        wcs=wcs, mask=mask)


def test_create_collection_from_spectrum1D():
    spec = Spectrum1D(spectral_axis=np.linspace(0, 50, 50) * u.AA,
                      flux=np.random.randn(50) * u.Jy,
                      uncertainty=StdDevUncertainty(np.random.sample(50), unit='Jy'))
    spec1 = Spectrum1D(spectral_axis=np.linspace(20, 60, 50) * u.AA,
                       flux=np.random.randn(50) * u.Jy,
                       uncertainty=StdDevUncertainty(np.random.sample(50), unit='Jy'))

    spec_coll = SpectrumCollection.from_spectra([spec, spec1])

    assert spec_coll.ndim == 1
    assert spec_coll.shape == (2, )
    assert spec_coll.nspectral == 50
    assert isinstance(spec_coll.flux, u.Quantity)
    assert isinstance(spec_coll.spectral_axis, u.Quantity)


def test_create_collection_from_spectra_without_uncertainties():
    spec = Spectrum1D(spectral_axis=np.linspace(0, 50, 50) * u.AA,
                      flux=np.random.randn(50) * u.Jy)
    spec1 = Spectrum1D(spectral_axis=np.linspace(20, 60, 50) * u.AA,
                       flux=np.random.randn(50) * u.Jy)

    spec_coll = SpectrumCollection.from_spectra([spec, spec1])