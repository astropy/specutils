import astropy.units as u
import numpy as np
import pytest
from astropy import time
from astropy.coordinates import (SkyCoord, EarthLocation, ICRS, Galactic,
                                 SpectralCoord, FK5)

from ..spectra.spectral_axis import SpectralAxis
from ..spectra.spectrum1d import Spectrum1D

from astropy.tests.helper import assert_quantity_allclose


def get_greenwich_earthlocation():
    """
    A helper function to get an EarthLocation for greenwich (without trying to
    do a download)
    """
    site_registry = EarthLocation._get_site_registry(force_builtin=True)
    return site_registry.get('greenwich')


def test_create_spectral_axis():

    site = get_greenwich_earthlocation()
    obstime = time.Time('2018-12-13 9:00')

    observer_gcrs = site.get_gcrs(obstime)

    wavelengths = np.linspace(500, 2500, 1001) * u.AA
    spectral_axis = SpectralAxis(wavelengths, observer=observer_gcrs)

    assert isinstance(spectral_axis, u.Quantity)
    assert len(spectral_axis) == 1001
    assert spectral_axis.bin_edges[0] == 499*u.AA


def test_create_with_bin_edges():

    wavelengths = np.linspace(500, 2500, 1001) * u.AA
    spectral_axis = SpectralAxis(wavelengths, bin_specification="edges")

    assert np.all(spectral_axis.bin_edges == wavelengths)
    assert spectral_axis[0] == 501*u.AA

    # Test irregular bin edges
    wavelengths = np.array([500, 510, 550, 560, 590])*u.AA
    spectral_axis = SpectralAxis(wavelengths, bin_specification="edges")

    assert np.all(spectral_axis.bin_edges == wavelengths)
    assert np.all(spectral_axis == [505., 530., 555., 575.]*u.AA)


# GENERAL TESTS

# We first run through a series of cases to test different ways of initializing
# the observer and target for SpectralAxis, including for example frames,
# SkyCoords, and making sure that SpectralAxis is not sensitive to the actual
# frame or representation class.

# Local Standard of Rest
LSRD = Galactic(u=0 * u.km, v=0 * u.km, w=0 * u.km,
                U=9 * u.km / u.s, V=12 * u.km / u.s, W=7 * u.km / u.s,
                representation_type='cartesian', differential_type='cartesian')

LSRD_EQUIV = [LSRD,
              SkyCoord(LSRD),  # as a SkyCoord
              LSRD.transform_to(ICRS()),  # different frame
              LSRD.transform_to(ICRS()).transform_to(Galactic())]  # different representation


@pytest.fixture(params=[None] + LSRD_EQUIV)
def observer(request):
    return request.param


# Target located in direction of motion of LSRD with no velocities
LSRD_DIR_STATIONARY = Galactic(u=9 * u.km, v=12 * u.km, w=7 * u.km,
                               representation_type='cartesian')

LSRD_DIR_STATIONARY_EQUIV = [
                             LSRD_DIR_STATIONARY,
                             SkyCoord(LSRD_DIR_STATIONARY),  # as a SkyCoord
                             LSRD_DIR_STATIONARY.transform_to(FK5()),  # different frame
                             LSRD_DIR_STATIONARY.transform_to(ICRS()).transform_to(Galactic())  # different representation
                            ]


@pytest.fixture(params=[None] + LSRD_DIR_STATIONARY_EQUIV)
def target(request):
    return request.param


def test_create_from_spectral_coord(observer, target):
    """
    Checks that parameters are correctly copied from the SpectralCoord object
    to the SpectralAxis object
    """
    spec_coord = SpectralCoord([100, 200, 300] * u.nm, observer=observer,
                               target=target, doppler_convention = 'optical',
                               doppler_rest = 6000*u.AA)
    spec_axis = SpectralAxis(spec_coord)
    assert spec_coord.observer == spec_axis.observer
    assert spec_coord.target == spec_axis.target
    assert spec_coord.radial_velocity == spec_axis.radial_velocity
    assert spec_coord.doppler_convention == spec_axis.doppler_convention
    assert spec_coord.doppler_rest == spec_axis.doppler_rest


def test_create_from_spectral_axis(observer, target):
    """
    Checks that parameters are correctly copied to the new SpectralAxis object
    """
    spec_axis1 = SpectralAxis([100, 200, 300] * u.nm, observer=observer,
                              target=target, doppler_convention = 'optical',
                              doppler_rest = 6000*u.AA)
    spec_axis2 = SpectralAxis(spec_axis1)
    assert spec_axis1.observer == spec_axis2.observer
    assert spec_axis1.target == spec_axis2.target
    assert spec_axis1.radial_velocity == spec_axis2.radial_velocity
    assert spec_axis1.doppler_convention == spec_axis2.doppler_convention
    assert spec_axis1.doppler_rest == spec_axis2.doppler_rest


def test_change_radial_velocity():
    wave = np.linspace(100, 200, 100) * u.AA
    flux = np.ones(100) * u.one
    spec = Spectrum1D(spectral_axis=wave, flux=flux,
                      radial_velocity=0 * u.km / u.s)

    assert spec.radial_velocity == 0 * u.km/u.s

    spec.set_radial_velocity_to(1 * u.km / u.s)

    assert spec.radial_velocity == 1 * u.km/u.s

    spec = Spectrum1D(spectral_axis=wave, flux=flux,
                      radial_velocity=10 * u.km / u.s)

    assert spec.radial_velocity == 10 * u.km / u.s

    spec.set_radial_velocity_to(5 * u.km / u.s)

    assert spec.radial_velocity == 5 * u.km / u.s


def test_no_change_radial_velocity():
    wave = np.linspace(100, 200, 100) * u.AA
    flux = np.ones(100) * u.one
    spec = Spectrum1D(spectral_axis=wave, flux=flux,
                      radial_velocity=0 * u.km / u.s)

    assert spec.radial_velocity == 0 * u.km/u.s
    spec.set_radial_velocity_to(10 * u.km/u.s)
    assert spec.radial_velocity == 10 * u.km/u.s
    assert_quantity_allclose(spec.wavelength, wave)


def test_change_redshift():
    wave = np.linspace(100, 200, 100) * u.AA
    flux = np.ones(100) * u.one
    spec = Spectrum1D(spectral_axis=wave, flux=flux, redshift=0)

    assert spec.redshift.unit.physical_type == 'dimensionless'
    assert_quantity_allclose(spec.redshift, u.Quantity(0))
    assert isinstance(spec.spectral_axis, SpectralAxis)

    spec.set_redshift_to(0.1)

    assert spec.redshift.unit.physical_type == 'dimensionless'
    assert_quantity_allclose(spec.redshift, u.Quantity(0.1))
    assert isinstance(spec.spectral_axis, SpectralAxis)

    spec = Spectrum1D(spectral_axis=wave, flux=flux, redshift=0.2)

    assert spec.redshift.unit.physical_type == 'dimensionless'
    assert_quantity_allclose(spec.redshift, u.Quantity(0.2))
    assert isinstance(spec.spectral_axis, SpectralAxis)

    spec.set_redshift_to(0.4)

    assert spec.redshift.unit.physical_type == 'dimensionless'
    assert_quantity_allclose(spec.redshift, u.Quantity(0.4))
    assert isinstance(spec.spectral_axis, SpectralAxis)


def test_no_change_redshift():
    wave = np.linspace(100, 200, 100) * u.AA
    flux = np.ones(100) * u.one
    spec = Spectrum1D(spectral_axis=wave, flux=flux, redshift=0)

    assert spec.redshift.unit.physical_type == 'dimensionless'
    assert_quantity_allclose(spec.redshift, u.Quantity(0))
    assert isinstance(spec.spectral_axis, SpectralAxis)

    spec.set_redshift_to(0.5)

    assert spec.redshift.unit.physical_type == 'dimensionless'
    assert_quantity_allclose(spec.redshift, u.Quantity(0.5))
    assert isinstance(spec.spectral_axis, SpectralAxis)

    assert_quantity_allclose(spec.wavelength, wave)


def test_pixel_descending_error():
    '''
    Spectral axes of pixel units must always be ascending.
    This test checks that an error is thrown if one is provided descending
    '''
    flux_unit = u.dimensionless_unscaled
    spec_unit = u.pix

    with pytest.raises(ValueError, match="u.pix spectral axes should always be ascending"):
        Spectrum1D(spectral_axis=(np.arange(5100, 5300)[::-1])*spec_unit,
                flux=np.random.randn(200)*flux_unit)
