import astropy.units as u
import numpy as np
import pytest
from astropy import time
from astropy.constants import c
from astropy.coordinates import (SkyCoord, EarthLocation, ICRS, GCRS, Galactic,
                                 CartesianDifferential,
                                 get_body_barycentric_posvel,
                                 FK5, CartesianRepresentation)

from ..extern.spectralcoord import SpectralCoord
from ..spectra.spectral_axis import SpectralAxis


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

LSRD_EQUIV = [
              LSRD,
              SkyCoord(LSRD),  # as a SkyCoord
              LSRD.transform_to(ICRS),  # different frame
              LSRD.transform_to(ICRS).transform_to(Galactic)  # different representation
              ]


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
