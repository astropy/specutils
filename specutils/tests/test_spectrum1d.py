import astropy.units as u
import astropy.wcs as fitswcs
import gwcs
import numpy as np
import pytest

from ..spectra.spectrum1d import Spectrum1D


def test_create_from_arrays():
    spec = Spectrum1D(spectral_axis=np.arange(50),
                      flux=np.random.randn(50))

    assert isinstance(spec.spectral_axis, u.Quantity)
    assert spec.spectral_axis.size == 50

    assert isinstance(spec.flux, u.Quantity)
    assert spec.flux.size == 50


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
    spec = Spectrum1D(spectral_axis=np.arange(50),
                      flux=np.random.randn(50))

    with pytest.raises(u.UnitsError) as e_info:
        spec.frequency

    spec = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                      flux=np.random.randn(49))

    assert spec.frequency.unit == u.GHz

    with pytest.raises(ValueError) as e_info:
        spec.velocity


def test_create_explicit_fitswcs():
    my_wcs = fitswcs.WCS(header={'CDELT1': 1, 'CRVAL1': 6562.8, 'CUNIT1': 'Angstrom',
                                 'CTYPE1': 'WAVE', 'RESTFRQ': 1400000000, 'CRPIX1': 25})

    spec = Spectrum1D(flux=[5,6,7] * u.Jy, wcs=my_wcs)
    spec.velocity_convention = "relativistic"

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
