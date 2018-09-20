import os

import astropy.units as u
import astropy.wcs as fitswcs
import gwcs
import numpy as np
import pytest
from astropy.nddata import StdDevUncertainty
from astropy.utils.data import get_pkg_data_filename
from astropy.modeling import models
from astropy.tests.helper import quantity_allclose
from ..manipulation import noise_region_uncertainty

from ..spectra import Spectrum1D, SpectralRegion


def test_create_from_arrays():
    spec = Spectrum1D(spectral_axis=np.arange(50) * u.AA,
                      flux=np.random.randn(50) * u.Jy)

    assert isinstance(spec.spectral_axis, u.Quantity)
    assert spec.spectral_axis.size == 50

    assert isinstance(spec.flux, u.Quantity)
    assert spec.flux.size == 50


def test_create_from_multidimensional_arrays():
    """
    This is a test for a bug that was fixed by #283. It makes sure that
    multidimensional flux arrays are handled properly when creating Spectrum1D
    objects.
    """

    freqs = np.arange(50) * u.GHz
    flux = np.random.random((5, len(freqs))) * u.Jy
    spec = Spectrum1D(spectral_axis=freqs, flux=flux)

    assert (spec.frequency == freqs).all()
    assert (spec.flux == flux).all()


def test_create_from_quantities():
    spec = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                      flux=np.random.randn(49) * u.Jy)

    assert isinstance(spec.spectral_axis, u.Quantity)
    assert spec.spectral_axis.unit == u.nm
    assert spec.spectral_axis.size == 49


def test_create_implicit_wcs():
    spec = Spectrum1D(spectral_axis=np.arange(50) * u.AA,
                      flux=np.random.randn(50) * u.Jy)

    assert isinstance(spec.wcs.wcs, gwcs.wcs.WCS)

    pix2world = spec.wcs.pixel_to_world(np.arange(5, 10))

    assert pix2world.size == 5
    assert isinstance(pix2world, np.ndarray)


def test_create_implicit_wcs_with_spectral_unit():
    spec = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                      flux=np.random.randn(49) * u.Jy)

    assert isinstance(spec.wcs.wcs, gwcs.wcs.WCS)

    pix2world = spec.wcs.pixel_to_world(np.arange(5, 10))

    assert pix2world.size == 5
    assert isinstance(pix2world, np.ndarray)


def test_spectral_axis_conversions():
    # By default the spectral axis units should be set to angstroms
    spec = Spectrum1D(flux=np.array([26.0, 44.5]) * u.Jy,
                      spectral_axis=np.array([400, 500]) * u.AA)

    assert np.all(spec.spectral_axis == np.array([400, 500]) * u.angstrom)
    assert spec.spectral_axis.unit == u.angstrom

    spec = Spectrum1D(spectral_axis=np.arange(50) * u.AA,
                      flux=np.random.randn(50) * u.Jy)

    assert spec.wavelength.unit == u.AA

    spec = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                      flux=np.random.randn(49) * u.Jy)

    assert spec.frequency.unit == u.GHz

    with pytest.raises(ValueError) as e_info:
        spec.velocity

    spec = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                      flux=np.random.randn(49) * u.Jy)

    new_spec = spec.with_spectral_unit(u.GHz)


def test_flux_unit_conversion():
    # By default the flux units should be set to Jy
    s = Spectrum1D(flux=np.array([26.0, 44.5]) * u.Jy,
                   spectral_axis=np.array([400, 500]) * u.nm)
    assert np.all(s.flux == np.array([26.0, 44.5]) * u.Jy)
    assert s.flux.unit == u.Jy

    # Simple Unit Conversion
    s = Spectrum1D(flux=np.array([26.0, 44.5]) * u.Jy, spectral_axis=np.array([400, 500])*u.nm)
    converted_spec = s.new_flux_unit(unit=u.uJy)[0]
    assert ((26.0 * u.Jy).to(u.uJy) == converted_spec.flux)

    # Make sure incompatible units raise UnitConversionError
    with pytest.raises(u.UnitConversionError):
        s.new_flux_unit(unit=u.m)

    # Pass custom equivalencies
    s = Spectrum1D(flux=np.array([26.0, 44.5]) * u.Jy,
                   spectral_axis=np.array([400, 500]) * u.nm)
    eq = [[u.Jy, u.m,
          lambda x: np.full_like(np.array(x), 1000.0, dtype=np.double),
          lambda x: np.full_like(np.array(x), 0.001, dtype=np.double)]]
    converted_spec = s.new_flux_unit(unit=u.m, equivalencies=eq)[0]
    assert 1000.0 * u.m == converted_spec.flux

    # Check if suppressing the unit conversion works
    s = Spectrum1D(flux=np.array([26.0, 44.5]) * u.Jy, spectral_axis=np.array([400, 500]) * u.nm)
    new_spec = s.new_flux_unit("uJy", suppress_conversion=True)
    assert new_spec.flux[0] == 26.0 * u.uJy


def test_wcs_transformations():
    # Test with a GWCS
    spec = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                      flux=np.random.randn(49) * u.Jy)

    pix_axis = spec.wcs.world_to_pixel(np.arange(20, 30))
    disp_axis = spec.wcs.pixel_to_world(np.arange(20, 30) * u.nm)

    assert isinstance(pix_axis, np.ndarray)
    assert isinstance(disp_axis, u.Quantity)

    # Test transform with different unit
    spec.wcs.world_to_pixel(np.arange(20, 30) * u.GHz)

    # Test with a FITS WCS
    my_wcs = fitswcs.WCS(header={'CDELT1': 1, 'CRVAL1': 6562.8, 'CUNIT1': 'Angstrom',
                                 'CTYPE1': 'WAVE', 'RESTFRQ': 1400000000, 'CRPIX1': 25})

    spec = Spectrum1D(flux=[5,6,7] * u.Jy, wcs=my_wcs)

    pix_axis = spec.wcs.world_to_pixel(np.arange(20, 30))
    disp_axis = spec.wcs.pixel_to_world(np.arange(20, 30) * u.nm)

    assert isinstance(pix_axis, np.ndarray)
    assert isinstance(disp_axis, u.Quantity)

def test_create_explicit_fitswcs():
    my_wcs = fitswcs.WCS(header={'CDELT1': 1, 'CRVAL1': 6562.8, 'CUNIT1': 'Angstrom',
                                 'CTYPE1': 'WAVE', 'RESTFRQ': 1400000000, 'CRPIX1': 25})

    spec = Spectrum1D(flux=[5,6,7] * u.Jy, wcs=my_wcs)
    spec = spec.with_velocity_convention("relativistic")

    assert isinstance(spec.spectral_axis, u.Quantity)
    assert spec.spectral_axis.unit.is_equivalent(u.AA)

    pix2world = spec.wcs.pixel_to_world(np.arange(3))

    assert pix2world.size == 3

    assert isinstance(spec.wavelength, u.Quantity)
    assert spec.wavelength.size == 3
    assert spec.wavelength.unit == u.AA

    assert isinstance(spec.frequency, u.Quantity)
    assert spec.frequency.size == 3
    assert spec.frequency.unit == u.GHz

    assert isinstance(spec.velocity, u.Quantity)
    assert spec.velocity.size == 3
    assert spec.velocity.unit == u.Unit('km/s')


def test_create_with_uncertainty():
    spec = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                      flux=np.random.sample(49) * u.Jy,
                      uncertainty=StdDevUncertainty(np.random.sample(49) * 0.1))

    assert isinstance(spec.uncertainty, StdDevUncertainty)

    spec = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                      flux=np.random.sample(49) * u.Jy,
                      uncertainty=StdDevUncertainty(np.random.sample(49) * 0.1))

    assert spec.flux.unit == spec.uncertainty.unit


def test_read_linear_solution():
    file_path = get_pkg_data_filename('data/L5g_0355+11_Cruz09.fits')

    spec = Spectrum1D.read(file_path, format='wcs1d-fits')

    assert isinstance(spec, Spectrum1D)

    assert isinstance(spec.flux, u.Quantity)
    assert isinstance(spec.spectral_axis, u.Quantity)

    assert spec.flux.size == spec.data.size
    assert spec.spectral_axis.size == spec.data.size


def test_energy_photon_flux():
    spec = Spectrum1D(spectral_axis=np.linspace(100, 1000, 10) * u.nm,
                      flux=np.random.randn(10)*u.Jy)
    assert spec.energy.size == 10
    assert spec.photon_flux.size == 10
    assert spec.photon_flux.unit == u.photon * u.cm**-2 * u.s**-1 * u.nm**-1


def test_repr():
    spec_with_wcs = Spectrum1D(spectral_axis=np.linspace(100, 1000, 10) * u.nm,
                               flux=np.random.random(10) * u.Jy)
    result = repr(spec_with_wcs)
    assert result.startswith('<Spectrum1D(flux=<Quantity [')
    assert 'spectral_axis=<Quantity [' in result

    spec_without_wcs = Spectrum1D(np.random.random(10) * u.Jy)
    result = repr(spec_without_wcs)
    assert result == '<Spectrum1D(flux={})>'.format(repr(spec_without_wcs.flux))


def test_str():
    spec = Spectrum1D(spectral_axis=np.linspace(100, 1000, 10) * u.nm,
                               flux=np.random.random(10) * u.Jy)
    result = str(spec)
    # Sanity check for contents of string representation
    assert result.startswith('Spectrum1D (length={})'.format(len(spec.flux)))
    lines = result.split('\n')
    flux = spec.flux
    sa = spec.spectral_axis
    assert len(lines) == 3
    assert lines[1].startswith('flux:')
    assert '[ {:.5}, ..., {:.5} ]'.format(flux[0], flux[-1]) in lines[1]
    assert 'mean={:.5}'.format(np.mean(flux)) in lines[1]
    assert lines[2].startswith('spectral axis:')
    assert '[ {:.5}, ..., {:.5} ]'.format(sa[0], sa[-1]) in lines[2]
    assert 'mean={:5}'.format(np.mean(sa)) in lines[2]

    # Test string representation with uncertainty
    spec_with_uncertainty = Spectrum1D(spectral_axis=np.linspace(100, 1000, 10) * u.nm,
                                       flux=np.random.random(10) * u.Jy,
                                       uncertainty=np.random.random(10))
    result = str(spec_with_uncertainty)
    lines = result.split('\n')
    unc = spec_with_uncertainty.uncertainty
    print(spec_with_uncertainty)
    assert len(lines) == 4
    assert lines[3].startswith('uncertainty')
    assert '[ {}, ..., {} ]'.format(unc[0], unc[-1]) in lines[3]

    # Test string representation with multiple flux
    spec_multi_flux = Spectrum1D(spectral_axis=np.linspace(100, 1000, 10) * u.nm,
                                 flux=np.random.random((3,10)) * u.Jy)
    result = str(spec_multi_flux)
    lines = result.split('\n')
    assert len(lines) == 5
    for i, line in enumerate(lines[1:4]):
        assert line.startswith('flux{:2}:'.format(i))

    # Test string representation with single-dimensional flux
    spec_single_flux = Spectrum1D(1 * u.Jy)
    result = str(spec_single_flux)
    assert result == 'Spectrum1D (length=1)\nflux:   {}'.format(spec_single_flux.flux)


def test_noise_estimate_uncertainty():

    np.random.seed(42)

    # Create values for spectrum.
    frequencies = np.linspace(1, 100, 10000) * u.um
    g = models.Gaussian1D(amplitude=1*u.Jy, mean=10*u.um, stddev=1*u.um)
    noise = np.random.normal(0., 0.01, frequencies.shape) * u.Jy
    flux = g(frequencies) + noise

    # Create the spectrum, spectral region and assign uncertainties
    spectrum = Spectrum1D(spectral_axis=frequencies, flux=flux)
    spectral_region = SpectralRegion(50*u.um, 80*u.um)
    spectrum_with_uncertainty = noise_region_uncertainty(spectrum, spectral_region)

    # Compute the expected.
    indices = np.nonzero((frequencies >= 50*u.um) & (frequencies <= 80*u.um))
    expected_uncertainty = np.std(flux[indices])*np.ones(len(frequencies))

    assert quantity_allclose(spectrum_with_uncertainty.uncertainty.array,
                             expected_uncertainty.value)

    # Same idea, but now with variance.  Create the spectrum, spectral region and assign uncertainties
    spectrum_with_uncertainty = noise_region_uncertainty(spectrum, spectral_region, np.var)

    # Compute the expected.
    indices = np.nonzero((frequencies >= 50*u.um) & (frequencies <= 80*u.um))
    expected_uncertainty = np.var(flux[indices])*np.ones(len(frequencies))

    assert quantity_allclose(spectrum_with_uncertainty.uncertainty.array,
                             expected_uncertainty.value)
