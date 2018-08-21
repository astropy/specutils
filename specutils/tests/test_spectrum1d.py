import astropy.units as u
import astropy.wcs as fitswcs
import gwcs
import numpy as np
import os
import pytest
from astropy.nddata import StdDevUncertainty
from astropy.utils.data import get_pkg_data_filename

from ..spectra.spectrum1d import Spectrum1D


def test_create_from_arrays():
    spec = Spectrum1D(spectral_axis=np.arange(50),
                      flux=np.random.randn(50))

    assert isinstance(spec.spectral_axis, u.Quantity)
    assert spec.spectral_axis.size == 50

    assert isinstance(spec.flux, u.Quantity)
    assert spec.flux.size == 50


def test_create_from_multidimensional_arrays():

    freqs = np.arange(50) * u.GHz
    flux = np.random.random((5, len(freqs))) * u.Jy
    spec = Spectrum1D(spectral_axis=freqs, flux=flux)

    assert (spec.frequency == freqs).all()
    assert (spec.flux == flux).all()


def test_create_from_quantities():
    spec = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                      flux=np.random.randn(49))

    assert isinstance(spec.spectral_axis, u.Quantity)
    assert spec.spectral_axis.unit == u.nm
    assert spec.spectral_axis.size == 49


def test_create_implicit_wcs():
    spec = Spectrum1D(spectral_axis=np.arange(50),
                      flux=np.random.randn(50))

    assert isinstance(spec.wcs.wcs, gwcs.wcs.WCS)

    pix2world = spec.wcs.pixel_to_world(np.arange(5, 10))

    assert pix2world.size == 5


def test_create_implicit_wcs_with_spectral_unit():
    spec = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                      flux=np.random.randn(49))

    assert isinstance(spec.wcs.wcs, gwcs.wcs.WCS)

    pix2world = spec.wcs.pixel_to_world(np.arange(5, 10))

    assert pix2world.size == 5


def test_spectral_axis_conversions():
    # By default the spectral axis units should be set to angstroms
    spec = Spectrum1D(flux=np.array([26.0, 44.5]) * u.Jy, spectral_axis=np.array([400, 500]))
    assert np.all(spec.spectral_axis == np.array([400, 500]) * u.angstrom)
    assert spec.spectral_axis.unit == u.angstrom

    spec = Spectrum1D(spectral_axis=np.arange(50),
                      flux=np.random.randn(50))

    assert spec.wavelength.unit == u.AA

    spec = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                      flux=np.random.randn(49))

    assert spec.frequency.unit == u.GHz

    with pytest.raises(ValueError) as e_info:
        spec.velocity

    spec = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                      flux=np.random.randn(49))

    new_spec = spec.with_spectral_unit(u.GHz)


def test_flux_unit_conversion():
    # By default the flux units should be set to Jy
    s = Spectrum1D(flux=np.array([26.0, 44.5]), spectral_axis=np.array([400, 500]) * u.nm)
    assert np.all(s.flux == np.array([26.0, 44.5]) * u.Jy)
    assert s.flux.unit == u.Jy

    # Simple Unit Conversion
    s = Spectrum1D(flux=np.array([26.0, 44.5]) * u.Jy, spectral_axis=np.array([400, 500])*u.nm)
    converted_value = s.to_flux(unit=u.uJy)[0]
    assert ((26.0 * u.Jy).to(u.uJy) == converted_value)

    # Make sure incompatible units raise UnitConversionError
    with pytest.raises(u.UnitConversionError):
        converted_value = s.to_flux(unit=u.m)

    # Pass custom equivalencies
    s = Spectrum1D(flux=np.array([26.0, 44.5]) * u.Jy, spectral_axis=np.array([400, 500]) * u.nm)
    eq = [[u.Jy, u.m,
          lambda x: np.full_like(np.array(x), 1000.0, dtype=np.double),
          lambda x: np.full_like(np.array(x), 0.001, dtype=np.double)]]
    converted_value = s.to_flux(unit=u.m, equivalencies=eq)[0]
    assert 1000.0 * u.m == converted_value

    # Check if suppressing the unit conversion works
    s = Spectrum1D(flux=np.array([26.0, 44.5]) * u.Jy, spectral_axis=np.array([400, 500]) * u.nm)
    s.to_flux("uJy", suppress_conversion=True)
    assert s.flux[0] == 26.0 * u.uJy


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
                      flux=np.random.sample(49),
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
