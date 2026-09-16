import warnings
from copy import deepcopy

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from astropy.tests.helper import assert_quantity_allclose

from ..spectra.spectral_axis import SpectralAxis, observer_for_frame
from ..spectra.spectral_frame import SpectralMedium, SPECTRAL_FRAMES, normalize_frame


def test_spectral_medium_vacuum():
    medium = SpectralMedium('vacuum')
    assert medium.is_vacuum and not medium.is_air
    assert medium.refraction_method is None
    assert medium.refraction_kwargs == {}
    assert str(medium) == 'vacuum'
    assert repr(medium) == "SpectralMedium(kind='vacuum')"
    assert medium == SpectralMedium('VACUUM')

    with pytest.raises(ValueError, match="Air conditions cannot be set"):
        SpectralMedium('vacuum', co2=400)


def test_spectral_medium_air():
    medium = SpectralMedium('air')
    assert medium.is_air
    assert medium.refraction_method == 'Morton2000'
    assert medium.refraction_kwargs['method'] == 'Morton2000'

    medium = SpectralMedium('air', refraction_method='Ciddor1996', co2=400,
                            temperature=10 * u.deg_C, pressure=700 * u.hPa,
                            humidity=20 * u.percent)
    assert medium.co2 == 400
    assert medium.temperature == 10 * u.deg_C
    assert medium.pressure == 700 * u.hPa
    assert medium.humidity == pytest.approx(0.2)
    assert medium.to_dict()['refraction_method'] == 'Ciddor1996'
    assert medium != SpectralMedium('air')


def test_spectral_medium_validation():
    with pytest.raises(ValueError, match="kind must be one of"):
        SpectralMedium('water')
    with pytest.raises(u.UnitsError):
        SpectralMedium('air', temperature=10 * u.m)
    with pytest.raises(u.UnitsError):
        SpectralMedium('air', pressure=10 * u.K)
    with pytest.raises(ValueError, match="humidity must be between"):
        SpectralMedium('air', humidity=1.5)


def test_spectral_medium_from_input():
    assert SpectralMedium.from_input(None) is None
    assert SpectralMedium.from_input('air') == SpectralMedium('air')
    medium = SpectralMedium('air', co2=400)
    assert SpectralMedium.from_input(medium) is medium
    assert SpectralMedium.from_input({'kind': 'air', 'co2': 400}) == medium
    with pytest.raises(TypeError):
        SpectralMedium.from_input(3)


@pytest.mark.parametrize('frame, expected', [
    (None, None),
    ('BARYCENT', 'BARYCENT'),
    ('barycent', 'BARYCENT'),
    ('barycentric', 'BARYCENT'),
    ('topocentric', 'TOPOCENT'),
    ('rest', 'SOURCE'),
    ('source', 'SOURCE'),
    ('lsrk', 'LSRK'),
])
def test_normalize_frame(frame, expected):
    assert normalize_frame(frame) == expected


def test_normalize_frame_invalid():
    with pytest.raises(ValueError, match="Unknown spectral frame"):
        normalize_frame('ecliptic')
    with pytest.raises(TypeError):
        normalize_frame(3)
    for frame in SPECTRAL_FRAMES:
        assert normalize_frame(frame) == frame


# ---------------------------------------------------------------------------
# SpectralAxis metadata
# ---------------------------------------------------------------------------


@pytest.fixture
def apo():
    return EarthLocation(lat=32.78 * u.deg, lon=-105.82 * u.deg, height=2788 * u.m)


@pytest.fixture
def obstime():
    return Time('2024-03-01T05:00:00')


@pytest.fixture
def target():
    return SkyCoord(ra=120 * u.deg, dec=-30 * u.deg, frame='icrs')


def test_spectral_axis_metadata(apo, obstime):
    axis = SpectralAxis(np.linspace(5000, 5010, 11) * u.AA, medium='air',
                        frame='barycentric', obstime='2024-03-01T05:00:00', location=apo)
    assert axis.medium == SpectralMedium('air')
    assert axis.frame == 'BARYCENT'
    assert axis.obstime == obstime
    assert axis.location is apo
    assert axis._metadata == {'medium': SpectralMedium('air'), 'frame': 'BARYCENT',
                              'obstime': obstime, 'location': apo}

    # Metadata survives slicing, copying, unit conversion and replication
    for other in (axis[2:5], axis.copy(), deepcopy(axis), axis.to(u.nm),
                  axis.replicate(value=np.arange(11.) * u.AA)):
        assert other.medium == axis.medium
        assert other.frame == axis.frame
        assert other.obstime == axis.obstime
        assert other.location is axis.location

    # ...and is inherited when wrapping an existing axis unless overridden
    wrapped = SpectralAxis(axis, medium='vacuum')
    assert wrapped.medium.is_vacuum
    assert wrapped.frame == 'BARYCENT'

    replaced = axis.replicate(frame='source', medium={'kind': 'air', 'co2': 400})
    assert replaced.frame == 'SOURCE'
    assert replaced.medium.co2 == 400


def test_spectral_axis_air_unit_guard():
    axis = SpectralAxis([5000., 5001.] * u.AA, medium='air')
    with pytest.raises(u.UnitConversionError, match="only holds in vacuum"):
        axis.to(u.GHz)
    with pytest.raises(u.UnitConversionError):
        axis.to(u.km / u.s, doppler_rest=5000 * u.AA, doppler_convention='optical')
    # wavelength to wavelength is fine, as is anything in vacuum
    assert axis.to(u.nm).unit == u.nm
    SpectralAxis([5000., 5001.] * u.AA, medium='vacuum').to(u.GHz)
    with pytest.raises(u.UnitsError, match="only be set on a wavelength axis"):
        SpectralAxis([1., 2.] * u.GHz, medium='air')
    with pytest.raises(u.UnitsError):
        SpectralAxis([1., 2.] * u.GHz, medium='vacuum').replicate(medium='air')
    # pixel axes are exempt
    SpectralAxis(np.arange(3) * u.pix, medium='air')


def test_spectral_axis_observer_from_location(apo, obstime, target):
    axis = SpectralAxis([5000.] * u.AA, frame='TOPOCENT', obstime=obstime, location=apo)
    assert axis.observer is not None
    assert axis.observer.obstime == obstime
    assert axis.target is None

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        topo = SpectralAxis([5000.] * u.AA, frame='TOPOCENT', obstime=obstime,
                            location=apo, target=target)
        bary = SpectralAxis([5000.] * u.AA, frame='BARYCENT', obstime=obstime,
                            location=apo, target=target)
        source = SpectralAxis([5000.] * u.AA, frame='SOURCE', obstime=obstime,
                              location=apo, target=target)
    # A target with no velocity: the topocentric radial velocity is entirely
    # the observer's motion, and vanishes in the barycentric and rest frames
    assert_quantity_allclose(topo.radial_velocity, 9.41 * u.km / u.s, atol=0.01 * u.km / u.s)
    assert_quantity_allclose(bary.radial_velocity, 0 * u.km / u.s, atol=1e-6 * u.km / u.s)
    assert_quantity_allclose(source.radial_velocity, 0 * u.km / u.s, atol=1e-6 * u.km / u.s)

    # No observer is built without a frame, or for frames astropy cannot represent
    assert SpectralAxis([5000.] * u.AA, obstime=obstime, location=apo).observer is None
    assert SpectralAxis([5000.] * u.AA, frame='CMBDIPOL', obstime=obstime,
                        location=apo).observer is None
    assert observer_for_frame('SOURCE', apo, obstime) is None

    with pytest.raises(ValueError, match="Set the velocity of the source on the target"):
        SpectralAxis([5000.] * u.AA, frame='TOPOCENT', obstime=obstime, location=apo,
                     target=target, radial_velocity=10 * u.km / u.s)


def test_spectral_axis_metadata_validation(apo):
    with pytest.raises(ValueError, match="Unknown spectral frame"):
        SpectralAxis([5000.] * u.AA, frame='nope')
    with pytest.raises(TypeError, match="EarthLocation"):
        SpectralAxis([5000.] * u.AA, location=(1, 2, 3))
    with pytest.raises(ValueError):
        SpectralAxis([5000.] * u.AA, obstime='not a time')
