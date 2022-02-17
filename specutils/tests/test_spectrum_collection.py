import astropy.units as u
import numpy as np
from numpy.testing import assert_allclose
import pytest
from astropy.nddata import StdDevUncertainty
from astropy.coordinates import SpectralCoord
from gwcs.wcs import WCS as GWCS

from ..spectra.spectrum1d import Spectrum1D
from ..spectra.spectrum_collection import SpectrumCollection
from ..utils.wcs_utils import gwcs_from_array


@pytest.fixture
def spectrum_collection():
    flux = u.Quantity(np.random.sample((5, 10)), unit='Jy')
    spectral_axis = u.Quantity(np.arange(50).reshape((5, 10)) + 1, unit='AA')
    wcs = np.array([gwcs_from_array(x) for x in spectral_axis])
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
    assert isinstance(single_spec.wcs, GWCS)
    assert single_spec.mask.shape == single_spec.flux.shape
    assert len(single_spec.meta) == len(spectrum_collection.meta[0])
    assert isinstance(single_spec.uncertainty, StdDevUncertainty)


def test_collection_without_optional_arguments():
    # Without uncertainties
    flux = u.Quantity(np.random.sample((5, 10)), unit='Jy')
    spectral_axis = u.Quantity(np.arange(50).reshape((5, 10)) + 1, unit='AA')
    uncertainty = StdDevUncertainty(np.random.sample((5, 10)), unit='Jy')
    wcs = np.array([gwcs_from_array(x) for x in spectral_axis])
    mask = np.ones((5, 10)).astype(bool)
    meta = [{'test': 5, 'info': [1, 2, 3]} for i in range(5)]

    spec_coll = SpectrumCollection(
        flux=flux, spectral_axis=spectral_axis, wcs=wcs, mask=mask, meta=meta)

    # Without mask
    spec_coll = SpectrumCollection(
        flux=flux, spectral_axis=spectral_axis, uncertainty=uncertainty,
        wcs=wcs, meta=meta)

    # Without meta
    spec_coll = SpectrumCollection(  # noqa
        flux=flux, spectral_axis=spectral_axis, uncertainty=uncertainty,
        wcs=wcs, mask=mask)


@pytest.mark.filterwarnings('ignore:Not all spectra have associated masks')
def test_create_collection_from_spectrum1D():
    spec = Spectrum1D(spectral_axis=SpectralCoord(np.linspace(0, 50, 50) * u.AA,
                      redshift=0.1), flux=np.random.randn(50) * u.Jy,
                      uncertainty=StdDevUncertainty(np.random.sample(50), unit='Jy'))
    spec1 = Spectrum1D(spectral_axis=SpectralCoord(np.linspace(20, 60, 50) * u.AA,
                       redshift=0.1), flux=np.random.randn(50) * u.Jy,
                       uncertainty=StdDevUncertainty(np.random.sample(50), unit='Jy'))

    spec_coll = SpectrumCollection.from_spectra([spec, spec1])

    assert spec_coll.ndim == 1
    assert spec_coll.shape == (2, )
    assert spec_coll.nspectral == 50
    assert isinstance(spec_coll.flux, u.Quantity)
    assert isinstance(spec_coll.spectral_axis, SpectralCoord)
    assert spec.spectral_axis.unit == spec_coll.spectral_axis.unit
    assert spec.flux.unit == spec_coll.flux.unit
    assert_allclose(spec_coll.spectral_axis.redshift.value, 0.1)


@pytest.mark.filterwarnings('ignore:Not all spectra have associated masks')
def test_create_collection_from_collections():
    spec = Spectrum1D(spectral_axis=np.linspace(0, 50, 50) * u.AA,
                      flux=np.random.randn(50) * u.Jy,
                      uncertainty=StdDevUncertainty(np.random.sample(50), unit='Jy'))
    spec1 = Spectrum1D(spectral_axis=np.linspace(20, 60, 50) * u.AA,
                       flux=np.random.randn(50) * u.Jy,
                       uncertainty=StdDevUncertainty(np.random.sample(50), unit='Jy'))

    spec_coll1 = SpectrumCollection.from_spectra([spec, spec1])

    spec2 = Spectrum1D(spectral_axis=np.linspace(40, 80, 50) * u.AA,
                       flux=np.random.randn(50) * u.Jy,
                       uncertainty=StdDevUncertainty(np.random.sample(50), unit='Jy'))

    spec_coll2 = SpectrumCollection.from_spectra([spec, spec2])

    spec_coll = SpectrumCollection.from_spectra([spec_coll1, spec_coll2, spec_coll1])

    assert spec_coll.ndim == 2
    assert spec_coll.shape == (3, 2)
    assert spec_coll.nspectral == 50
    assert isinstance(spec_coll.flux, u.Quantity)
    assert isinstance(spec_coll.spectral_axis, u.Quantity)
    assert spec.spectral_axis.unit == spec_coll.spectral_axis.unit
    assert spec.flux.unit == spec_coll.flux.unit


@pytest.mark.filterwarnings('ignore:Not all spectra have associated uncertainties of the same type')
@pytest.mark.filterwarnings('ignore:Not all spectra have associated masks')
def test_create_collection_from_spectra_without_uncertainties():
    spec = Spectrum1D(spectral_axis=np.linspace(0, 50, 50) * u.AA,
                      flux=np.random.randn(50) * u.Jy)
    spec1 = Spectrum1D(spectral_axis=np.linspace(20, 60, 50) * u.AA,
                       flux=np.random.randn(50) * u.Jy)

    SpectrumCollection.from_spectra([spec, spec1])


def test_mismatched_spectral_axes_parameters():
    spec = Spectrum1D(spectral_axis=SpectralCoord(np.linspace(0, 50, 50) * u.AA,
                      radial_velocity=u.Quantity(100.0, "km/s")),
                      flux=np.random.randn(50) * u.Jy)
    spec1 = Spectrum1D(spectral_axis=SpectralCoord(np.linspace(20, 60, 50) * u.AA,
                       radial_velocity=u.Quantity(200.0, "km/s")),
                       flux=np.random.randn(50) * u.Jy)

    with pytest.raises(ValueError):
        SpectrumCollection.from_spectra([spec, spec1])


@pytest.mark.parametrize('scshape,expected_len', [((5, 10), 5), ((4, 5, 10), 4)])
def test_len(scshape, expected_len):
    flux = u.Quantity(np.random.sample(scshape), unit='Jy')
    spectral_axis = u.Quantity(np.arange(np.prod(scshape)).reshape(scshape) + 1, unit='AA')
    sc2d = SpectrumCollection(flux=flux, spectral_axis=spectral_axis)

    assert sc2d.shape == scshape[:-1]
    assert len(sc2d) == expected_len
