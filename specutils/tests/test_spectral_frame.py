import warnings
from copy import deepcopy

import astropy.units as u
import numpy as np
import pytest
from astropy.constants import c
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.nddata import StdDevUncertainty
from astropy.time import Time
from astropy.utils.exceptions import AstropyUserWarning
from astropy.tests.helper import assert_quantity_allclose

from ..spectra.spectral_axis import SpectralAxis, observer_for_frame
from ..spectra.spectrum import Spectrum
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


# ---------------------------------------------------------------------------
# Spectrum metadata
# ---------------------------------------------------------------------------

@pytest.fixture
def spectrum(apo, obstime, target):
    wavelength = np.linspace(5000, 5010, 11) * u.AA
    flux = np.ones(11) * u.Jy
    uncertainty = StdDevUncertainty(0.1 * np.ones(11) * u.Jy)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return Spectrum(spectral_axis=wavelength, flux=flux, uncertainty=uncertainty,
                        medium='air', frame='TOPOCENT', obstime=obstime, location=apo,
                        target=target)


def test_spectrum_metadata(spectrum, apo, obstime):
    assert spectrum.medium == SpectralMedium('air')
    assert spectrum.frame == 'TOPOCENT'
    assert spectrum.in_rest_frame is False
    assert spectrum.obstime == obstime
    assert spectrum.location is apo
    assert spectrum.observer is not None
    assert spectrum.target is not None
    # Topocentric radial velocity of a target with no velocity is the
    # observer's own motion
    assert_quantity_allclose(spectrum.radial_velocity, 9.41 * u.km / u.s,
                             atol=0.01 * u.km / u.s)
    assert "medium=air; frame=TOPOCENT" in repr(spectrum)
    assert "Medium=air\nFrame=TOPOCENT" in str(spectrum)

    plain = Spectrum(spectral_axis=spectrum.spectral_axis.quantity, flux=spectrum.flux)
    assert plain.medium is None and plain.frame is None and plain.in_rest_frame is None
    assert plain.target is None and plain.observer is None
    assert plain.obstime is None and plain.location is None
    assert "medium" not in repr(plain)


def test_spectrum_metadata_propagates(spectrum):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        derived = {
            'slice': spectrum[2:5],
            'world slice': spectrum[5002 * u.AA:5006 * u.AA],
            'multiply': spectrum * 2,
            'add': spectrum + 1 * u.Jy,
            'subtract': spectrum - spectrum,
            'divide': spectrum / spectrum,
            'power': spectrum ** 2,
            'copy': spectrum._copy(),
            'spectral unit': spectrum.with_spectral_axis_unit(u.nm),
            'flux unit': spectrum.with_flux_unit(u.mJy),
            'velocity convention': spectrum.with_velocity_convention('optical'),
        }
    for name, other in derived.items():
        assert other.medium == spectrum.medium, name
        assert other.frame == spectrum.frame, name
        assert other.obstime == spectrum.obstime, name
        assert other.location is spectrum.location, name
        assert other.observer is not None, name
        assert other.target is not None, name
        assert_quantity_allclose(other.radial_velocity, spectrum.radial_velocity)


def test_spectrum_metadata_from_spectral_axis(spectrum):
    axis = SpectralAxis(np.linspace(5000, 5010, 11) * u.AA, medium='air')
    flux = np.ones(11) * u.Jy

    # Metadata missing from the axis can be supplied alongside it...
    spec = Spectrum(spectral_axis=axis, flux=flux, frame='BARYCENT')
    assert spec.medium.is_air
    assert spec.frame == 'BARYCENT'

    # ...but metadata already on the axis cannot be overridden
    with pytest.raises(ValueError, match="Cannot separately set medium"):
        Spectrum(spectral_axis=axis, flux=flux, medium='vacuum')
    with pytest.raises(ValueError, match="Cannot separately set frame, target"):
        Spectrum(spectral_axis=spectrum.spectral_axis, flux=flux, frame='SOURCE',
                 target=spectrum.target)

    # The air guard applies to the Spectrum conveniences too
    with pytest.raises(u.UnitConversionError):
        spec.frequency
    with pytest.raises(u.UnitConversionError):
        spec.energy
    with pytest.raises(u.UnitConversionError):
        spec.with_spectral_axis_unit(u.GHz)
    assert spec.wavelength.unit == u.AA


def test_spectrum_metadata_cube():
    from astropy.wcs import WCS

    flux = np.arange(24).reshape([2, 3, 4]) * u.Jy
    wcs = WCS({"CTYPE1": "RA---TAN", "CTYPE2": "DEC--TAN", "CTYPE3": "WAVE-LOG",
               "CRVAL1": 205, "CRVAL2": 27, "CRVAL3": 3.622e-7,
               "CDELT1": -0.0001, "CDELT2": 0.0001, "CDELT3": 8e-11,
               "CRPIX1": 0, "CRPIX2": 0, "CRPIX3": 0})
    spec = Spectrum(flux=flux, wcs=wcs, frame='BARYCENT', medium='vacuum', redshift=0.01)
    assert spec.frame == 'BARYCENT'
    assert spec.medium.is_vacuum

    for other in (spec.with_spectral_axis_last(), spec[:, 1:, :], spec[:, 1, 2],
                  spec.mean(axis='spatial')):
        assert other.frame == 'BARYCENT'
        assert other.medium.is_vacuum
        assert_quantity_allclose(other.redshift, 0.01)


# ---------------------------------------------------------------------------
# Manipulation functions
# ---------------------------------------------------------------------------

def test_manipulation_keeps_metadata(spectrum):
    from ..manipulation import (FluxConservingResampler, LinearInterpolatedResampler,
                                SplineInterpolatedResampler, extract_region, excise_regions,
                                gaussian_smooth)
    from ..spectra.spectral_region import SpectralRegion

    grid = np.linspace(5001, 5009, 5) * u.AA
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        derived = {
            'flux conserving': FluxConservingResampler()(spectrum, grid),
            'linear': LinearInterpolatedResampler()(spectrum, grid),
            'spline': SplineInterpolatedResampler()(spectrum, grid),
            'truncated': LinearInterpolatedResampler('truncate')(
                spectrum, np.linspace(4990, 5009, 5) * u.AA),
            'extract': extract_region(spectrum, SpectralRegion(5002 * u.AA, 5006 * u.AA)),
            'extract joined': extract_region(
                spectrum, SpectralRegion([(5001 * u.AA, 5003 * u.AA),
                                          (5006 * u.AA, 5008 * u.AA)]),
                return_single_spectrum=True),
            'extract empty': extract_region(spectrum, SpectralRegion(6000 * u.AA, 6010 * u.AA)),
            'excise': excise_regions(spectrum, [SpectralRegion(5002 * u.AA, 5006 * u.AA)]),
            'smooth': gaussian_smooth(spectrum, 1),
        }
    for name, other in derived.items():
        assert other.medium == spectrum.medium, name
        assert other.frame == spectrum.frame, name
        assert other.obstime == spectrum.obstime, name
        assert other.observer is not None, name
        assert other.target is not None, name

    # A new grid carrying different metadata is interpreted in the frame of
    # the spectrum being resampled, with a warning
    grid = SpectralAxis(grid, frame='BARYCENT')
    with pytest.warns(AstropyUserWarning, match="interpreted in the frame of the input"):
        resampled = LinearInterpolatedResampler()(spectrum, grid)
    assert resampled.frame == 'TOPOCENT'


# ---------------------------------------------------------------------------
# Frame and medium conversions
# ---------------------------------------------------------------------------

KMS = u.km / u.s


@pytest.fixture
def moving_target():
    return SkyCoord(ra=120 * u.deg, dec=-30 * u.deg, radial_velocity=50 * KMS,
                    distance=100 * u.pc)


@pytest.fixture
def topocentric(apo, obstime, moving_target):
    wavelength = np.linspace(5000, 5010, 11) * u.AA
    return Spectrum(spectral_axis=wavelength, flux=np.ones(11) * u.Jy, medium='vacuum',
                    frame='TOPOCENT', obstime=obstime, location=apo, target=moving_target)


def _doppler(wavelength, velocity):
    beta = (velocity / c).to_value(u.one)
    return wavelength * np.sqrt((1 + beta) / (1 - beta))


def test_with_frame_from_observer_and_target(topocentric):
    spec = topocentric
    # Topocentric radial velocity is the source velocity plus the observer's motion
    bc = spec.barycentric_correction
    assert_quantity_allclose(bc, -9.41 * KMS, atol=0.01 * KMS)
    assert_quantity_allclose(spec.radial_velocity, 50 * KMS - bc)
    # ...which agrees with astropy's independent calculation to a few m/s
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', AstropyUserWarning)
        astropy_bc = SkyCoord(spec.target).radial_velocity_correction(
            'barycentric', obstime=spec.obstime, location=spec.location)
    assert_quantity_allclose(bc, astropy_bc, atol=0.01 * KMS)

    bary = spec.with_frame('BARYCENT')
    assert bary.frame == 'BARYCENT'
    assert_quantity_allclose(bary.radial_velocity, 50 * KMS, atol=1e-6 * KMS)
    assert_quantity_allclose(bary.barycentric_correction, 0 * KMS, atol=1e-6 * KMS)
    assert_quantity_allclose(bary.spectral_axis, _doppler(spec.spectral_axis, bc))
    assert 'original_wcs' in bary.meta

    rest = spec.to_rest()
    assert rest.frame == 'SOURCE'
    assert rest.in_rest_frame
    assert_quantity_allclose(rest.radial_velocity, 0 * KMS, atol=1e-6 * KMS)
    assert_quantity_allclose(rest.spectral_axis,
                             _doppler(spec.spectral_axis, -spec.radial_velocity))
    assert_quantity_allclose(bary.to_rest().spectral_axis, rest.spectral_axis)

    # Round trip back to the telescope frame
    back = rest.with_frame('TOPOCENT')
    assert back.frame == 'TOPOCENT'
    assert_quantity_allclose(back.spectral_axis, spec.spectral_axis)
    assert_quantity_allclose(back.radial_velocity, spec.radial_velocity)

    helio = spec.with_frame('heliocentric')
    assert helio.frame == 'HELIOCEN'
    assert_quantity_allclose(helio.radial_velocity, 50 * KMS, atol=0.02 * KMS)
    assert spec.with_frame('LSRK').frame == 'LSRK'

    same = spec.with_frame('TOPOCENT')
    assert same is not spec and same.frame == 'TOPOCENT'
    assert_quantity_allclose(same.spectral_axis, spec.spectral_axis)

    with pytest.raises(ValueError, match="astropy has no coordinate frame"):
        spec.with_frame('CMBDIPOL')
    with pytest.raises(ValueError, match="already in the 'TOPOCENT' frame"):
        spec.with_frame('TOPOCENT', velocity=1 * KMS)
    with pytest.raises(ValueError, match="A frame must be given"):
        spec.with_frame(None)


def test_with_frame_from_velocity(topocentric):
    reference = topocentric
    bc = reference.barycentric_correction
    spec = Spectrum(spectral_axis=reference.spectral_axis.quantity, flux=reference.flux,
                    frame='TOPOCENT', radial_velocity=reference.radial_velocity)
    assert spec.barycentric_correction is None

    bary = spec.with_frame('BARYCENT', velocity=bc)
    assert bary.frame == 'BARYCENT'
    assert_quantity_allclose(bary.radial_velocity, 50 * KMS, atol=1e-6 * KMS)
    assert_quantity_allclose(bary.spectral_axis, reference.with_frame('BARYCENT').spectral_axis)

    rest = bary.to_rest()
    assert rest.frame == 'SOURCE'
    assert_quantity_allclose(rest.radial_velocity, 0 * KMS, atol=1e-6 * KMS)
    assert_quantity_allclose(rest.spectral_axis, reference.to_rest().spectral_axis)

    with pytest.raises(ValueError, match="without either the ``velocity``"):
        spec.with_frame('BARYCENT')
    with pytest.raises(u.UnitsError):
        spec.with_frame('BARYCENT', velocity=3 * u.AA)
    with pytest.raises(ValueError, match="current frame of the spectrum is unknown"):
        Spectrum(spectral_axis=spec.spectral_axis.quantity, flux=spec.flux).to_rest()

    # Applying a velocity with an observer and target moves the observer
    mixed = reference.with_frame('BARYCENT', velocity=bc)
    assert_quantity_allclose(mixed.radial_velocity, 50 * KMS, atol=1e-6 * KMS)
    assert_quantity_allclose(mixed.barycentric_correction, 0 * KMS, atol=1e-6 * KMS)


def test_with_frame_rebuilds_observer(apo, obstime):
    wavelength = np.linspace(5000, 5010, 11) * u.AA
    spec = Spectrum(spectral_axis=wavelength, flux=np.ones(11) * u.Jy, frame='TOPOCENT',
                    obstime=obstime, location=apo, radial_velocity=10 * KMS)
    assert spec.observer is not None and spec.target is None

    bary = spec.with_frame('BARYCENT', velocity=-9.41 * KMS)
    assert bary.frame == 'BARYCENT'
    assert bary.observer.__class__.__name__ == 'ICRS'
    assert_quantity_allclose(bary.radial_velocity, 0.59 * KMS)

    # Without a location the observer cannot follow the frame and is dropped
    spec = Spectrum(spectral_axis=wavelength, flux=np.ones(11) * u.Jy, frame='TOPOCENT',
                    observer=apo.get_gcrs(obstime), radial_velocity=10 * KMS)
    with pytest.warns(AstropyUserWarning, match="Dropping the observer"):
        bary = spec.with_frame('BARYCENT', velocity=-9.41 * KMS)
    assert bary.observer is None


def test_with_medium(topocentric):
    from ..utils.wcs_utils import vac_to_air

    spec = topocentric
    air = spec.with_medium('air')
    assert air.medium == SpectralMedium('air')
    assert_quantity_allclose(air.spectral_axis, vac_to_air(spec.spectral_axis.quantity))
    # Frame and observer metadata are untouched
    assert air.frame == spec.frame
    assert_quantity_allclose(air.radial_velocity, spec.radial_velocity)
    assert 'original_wcs' in air.meta

    vacuum = air.with_medium('vacuum', scheme='iteration')
    assert vacuum.medium.is_vacuum
    assert_quantity_allclose(vacuum.spectral_axis, spec.spectral_axis, atol=1e-9 * u.AA)

    conditions = SpectralMedium('air', pressure=700 * u.hPa, temperature=5 * u.deg_C)
    thin_air = air.with_medium(conditions)
    assert thin_air.medium == conditions
    assert_quantity_allclose(thin_air.spectral_axis,
                             vac_to_air(spec.spectral_axis.quantity, pressure=700 * u.hPa,
                                        temperature=5 * u.deg_C))

    same = air.with_medium('air')
    assert same is not air
    assert_quantity_allclose(same.spectral_axis, air.spectral_axis)

    with pytest.warns(AstropyUserWarning, match="velocity shift to air wavelengths"):
        air.to_rest()

    with pytest.raises(ValueError, match="medium of this spectrum is unknown"):
        Spectrum(spectral_axis=spec.spectral_axis.quantity, flux=spec.flux).with_medium('air')
    with pytest.raises(u.UnitsError, match="must be in wavelength units"):
        spec.with_spectral_axis_unit(u.GHz).with_medium('air')
    with pytest.raises(ValueError, match="A medium must be given"):
        spec.with_medium(None)


# ---------------------------------------------------------------------------
# FITS keywords
# ---------------------------------------------------------------------------

@pytest.fixture
def linear_wcs_header():
    return {'CTYPE1': 'WAVE', 'CUNIT1': 'Angstrom', 'CRPIX1': 1, 'CRVAL1': 5000, 'CDELT1': 1}


def test_metadata_from_fits_wcs(linear_wcs_header):
    from astropy.wcs import WCS

    flux = np.arange(1, 11) * u.Jy
    header = dict(linear_wcs_header, SPECSYS='BARYCENT', **{'MJD-OBS': 60000.0, 'MJD-BEG': 60000.0,
                  'MJD-END': 60000.02, 'ZSOURCE': 0.001, 'TIMESYS': 'TAI',
                  'OBSGEO-X': -1463969.3, 'OBSGEO-Y': -5166673.3, 'OBSGEO-Z': 3434985.7})
    spec = Spectrum(flux=flux, wcs=WCS(header))
    assert spec.medium.is_vacuum
    assert spec.frame == 'BARYCENT'
    assert spec.obstime.scale == 'tai'
    assert_quantity_allclose(spec.obstime.mjd, 60000.01)
    assert_quantity_allclose(spec.location.geodetic.lat, 32.78 * u.deg, atol=1e-3 * u.deg)
    assert_quantity_allclose(spec.redshift, 0.001)

    # Air wavelengths, MJD-AVG preferred, ZSOURCE meaningless in the rest frame
    header = dict(linear_wcs_header, CTYPE1='AWAV', SPECSYS='SOURCE', ZSOURCE=0.001,
                  **{'MJD-AVG': 60000.5, 'MJD-OBS': 60000.0})
    spec = Spectrum(flux=flux, wcs=WCS(header))
    assert spec.medium.is_air
    assert spec.frame == 'SOURCE'
    assert spec.obstime.scale == 'utc' and spec.obstime.mjd == 60000.5
    assert spec.redshift == 0
    assert spec.location is None

    # Explicit arguments take precedence over the WCS
    spec = Spectrum(flux=flux, wcs=WCS(header), medium='vacuum', frame='TOPOCENT')
    assert spec.medium.is_vacuum and spec.frame == 'TOPOCENT'

    # Unknown SPECSYS is ignored with a warning; frequency axes are vacuum
    header = dict(linear_wcs_header, CTYPE1='FREQ', CUNIT1='GHz', SPECSYS='NOPE')
    with pytest.warns(AstropyUserWarning, match="Ignoring unrecognised SPECSYS"):
        spec = Spectrum(flux=flux, wcs=WCS(header))
    assert spec.frame is None and spec.medium.is_vacuum

    # Nothing is invented for a bare WCS
    spec = Spectrum(flux=flux, wcs=WCS(dict(linear_wcs_header, CTYPE1='VRAD', CUNIT1='km/s')))
    assert spec.medium is None and spec.frame is None and spec.obstime is None


@pytest.mark.parametrize('fmt', ['wcs1d-fits', 'tabular-fits'])
def test_fits_metadata_round_trip(tmp_path, fmt, linear_wcs_header, apo, obstime, target):
    from astropy.io import fits
    from astropy.wcs import WCS

    flux = np.arange(1, 11) * u.Jy
    spec = Spectrum(flux=flux, wcs=WCS(linear_wcs_header), medium='air', frame='TOPOCENT',
                    obstime=obstime, location=apo, target=target)
    path = tmp_path / 'spec.fits'
    if fmt == 'wcs1d-fits':
        spec.write(path, format=fmt, hdu=0)
    else:
        spec.write(path, format=fmt)

    with fits.open(path) as hdulist:
        header = fits.Header(hdulist[0].header)
        header.update(hdulist[-1].header)
    assert header['SPECSYS'] == 'TOPOCENT'
    assert header['TIMESYS'] == 'UTC'
    assert_quantity_allclose(header['MJD-AVG'], obstime.mjd)
    assert header['DATE-AVG'] == obstime.isot
    assert header['RADESYS'] == 'ICRS'
    assert_quantity_allclose(header['RA'], 120)
    assert_quantity_allclose(header['DEC'], -30)
    assert header['CTYPE1' if fmt == 'wcs1d-fits' else 'TCTYP1'] == 'AWAV'
    assert_quantity_allclose(header['OBSGEO-X'], apo.geocentric[0].to_value(u.m), rtol=1e-6)

    other = Spectrum.read(path, format=fmt)
    assert other.medium == spec.medium
    assert other.frame == spec.frame
    assert other.obstime.isot == obstime.isot
    assert_quantity_allclose(other.location.geodetic.height, apo.geodetic.height, atol=1 * u.m)
    assert_quantity_allclose(SkyCoord(other.target).separation(target), 0 * u.deg,
                             atol=1e-6 * u.deg)
    assert_quantity_allclose(other.radial_velocity, spec.radial_velocity, atol=1e-3 * KMS)
    assert_quantity_allclose(other.spectral_axis, spec.spectral_axis)

    # Spectra without metadata write none
    Spectrum(flux=flux, wcs=WCS(linear_wcs_header)).write(path, format=fmt, overwrite=True,
                                                          **({'hdu': 0} if 'wcs' in fmt else {}))
    with fits.open(path) as hdulist:
        for hdu in hdulist:
            assert not any(key in hdu.header for key in ('SPECSYS', 'MJD-AVG', 'RA', 'OBSGEO-X'))


def test_spectral_axis_metadata_from_header():
    from astropy.io import fits
    from ..io.parsing_utils import spectral_axis_metadata_from_header

    def read(**keywords):
        header = fits.Header()
        header.update(keywords)
        return spectral_axis_metadata_from_header(header)

    assert read() == {}

    # Observation time: the mid-point of the exposure
    assert read(**{'DATE-OBS': '2024-03-01T05:00:00', 'EXPTIME': 600.0})['obstime'].isot == \
        '2024-03-01T05:05:00.000'
    assert read(**{'DATE-OBS': '2024-03-01', 'TIME-OBS': '05:00:00',
                   'EXPTIME': 600.0})['obstime'].isot == '2024-03-01T05:05:00.000'
    assert read(**{'DATE-BEG': '2024-03-01T05:00:00',
                   'DATE-END': '2024-03-01T05:10:00'})['obstime'].isot == '2024-03-01T05:05:00.000'
    assert read(**{'MJD-BEG': 60000.0, 'MJD-END': 60000.02})['obstime'].mjd == 60000.01
    assert read(**{'DATE-OBS': '2024-03-01T05:00:00', 'EXPTIME': 600.0,
                   'MJD-AVG': 60000.5})['obstime'].mjd == 60000.5
    assert read(**{'MJD': 60000.0})['obstime'].mjd == 60000.0  # no EXPTIME: the start
    tai = read(**{'MJD-AVG': 60000.5, 'TIMESYS': 'TAI'})['obstime']
    assert tai.scale == 'tai' and tai.mjd == 60000.5
    assert 'obstime' not in read(**{'DATE-OBS': '01/03/24'})

    header = fits.Header()
    header['MIDTIME'] = 60000.25
    assert spectral_axis_metadata_from_header(header, time_key='MIDTIME')['obstime'].mjd == \
        60000.25

    # Target position in degrees or sexagesimal, in the RADESYS frame
    target = read(RA=120.0, DEC=-30.0)['target']
    assert target.frame.name == 'icrs'
    assert_quantity_allclose([target.ra.deg, target.dec.deg], [120, -30])
    target = read(RA='08:00:00.0', DEC='-30:00:00', RADESYS='FK5', EQUINOX=2000.0)['target']
    assert target.frame.name == 'fk5'
    assert_quantity_allclose([target.ra.deg, target.dec.deg], [120, -30])
    target = read(RA='08 00 00', DEC='-30 00 00')['target']
    assert_quantity_allclose([target.ra.deg, target.dec.deg], [120, -30])
    target = read(RA_TARG=121.0, DEC_TARG=-31.0, RA=1.0, DEC=1.0)['target']
    assert_quantity_allclose(target.ra.deg, 121)
    assert 'target' not in read(RA='N/A', DEC='N/A')
    header = fits.Header()
    header.update(dict(RA=1.0, DEC=1.0, MYRA=120.0, MYDEC=-30.0))
    target = spectral_axis_metadata_from_header(header, target_keys=('MYRA', 'MYDEC'))['target']
    assert_quantity_allclose(target.ra.deg, 120)

    # Medium and frame from the WCS keywords, or given explicitly
    assert read(NAXIS=1, CTYPE1='AWAV', SPECSYS='TOPOCENT') == {
        'medium': SpectralMedium('air'), 'frame': 'TOPOCENT'}
    assert read(NAXIS=2, CTYPE1='RA---TAN', CTYPE2='FREQ')['medium'].is_vacuum
    assert read(TCTYP1='WAVE')['medium'].is_vacuum
    assert 'medium' not in read(NAXIS=1, CTYPE1='VRAD')
    with pytest.warns(AstropyUserWarning, match="Ignoring unrecognised SPECSYS"):
        assert 'frame' not in read(SPECSYS='NOPE')
    header = fits.Header()
    header.update(dict(CTYPE1='AWAV', SPECSYS='TOPOCENT'))
    assert spectral_axis_metadata_from_header(header, medium='vacuum', frame='barycentric') == {
        'medium': SpectralMedium('vacuum'), 'frame': 'BARYCENT'}

    # Location from OBSGEO keywords or given explicitly
    location = read(**{'OBSGEO-X': -1463969.3, 'OBSGEO-Y': -5166673.3,
                       'OBSGEO-Z': 3434985.7})['location']
    assert_quantity_allclose(location.geodetic.lat, 32.78 * u.deg, atol=1e-3 * u.deg)
    location = read(**{'OBSGEO-L': -105.82, 'OBSGEO-B': 32.78, 'OBSGEO-H': 2788.0})['location']
    assert_quantity_allclose(location.geodetic.lon, -105.82 * u.deg)
    apo = EarthLocation(lat=32.78 * u.deg, lon=-105.82 * u.deg, height=2788 * u.m)
    assert spectral_axis_metadata_from_header(fits.Header(), location=apo)['location'] is apo
    with pytest.raises(TypeError):
        spectral_axis_metadata_from_header(fits.Header(), location=3)


# ---------------------------------------------------------------------------
# SpectrumCollection
# ---------------------------------------------------------------------------

def test_spectrum_collection_metadata(apo, obstime, moving_target):
    from ..spectra.spectrum_collection import SpectrumCollection

    flux = np.ones((3, 11)) * u.Jy
    spectral_axis = np.tile(np.linspace(5000, 5010, 11), (3, 1)) * u.AA
    collection = SpectrumCollection(flux, spectral_axis=spectral_axis, medium='vacuum',
                                    frame='TOPOCENT', obstime=obstime, location=apo,
                                    target=moving_target)
    assert collection.medium.is_vacuum
    assert collection.frame == 'TOPOCENT'
    assert collection.obstime == obstime
    assert collection.location is apo
    assert collection.observer is not None
    assert collection.in_rest_frame is False
    assert_quantity_allclose(collection.barycentric_correction, -9.41 * KMS, atol=0.01 * KMS)
    assert "Medium:              vacuum" in repr(collection)

    # Individual spectra carry the metadata, and can be recombined
    spec = collection[1]
    assert spec.frame == 'TOPOCENT' and spec.medium.is_vacuum and spec.target is not None
    assert_quantity_allclose(spec.radial_velocity, collection.radial_velocity)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)  # no uncertainties or masks
        rebuilt = SpectrumCollection.from_spectra([collection[i] for i in range(3)])
    assert rebuilt.frame == 'TOPOCENT' and rebuilt.medium.is_vacuum
    assert rebuilt.obstime == obstime and rebuilt.location is apo
    assert_quantity_allclose(rebuilt.radial_velocity, collection.radial_velocity)

    # Frame and medium conversions work on the whole collection
    rest = collection.to_rest()
    assert rest.frame == 'SOURCE'
    assert_quantity_allclose(rest.radial_velocity, 0 * KMS, atol=1e-6 * KMS)
    assert_quantity_allclose(rest.spectral_axis[1], collection[1].to_rest().spectral_axis)
    air = collection.with_medium('air')
    assert air.medium.is_air
    assert_quantity_allclose(air.spectral_axis[0], collection[0].with_medium('air').spectral_axis)

    # Mixed metadata cannot be combined
    topo = Spectrum(spectral_axis=spectral_axis[0], flux=flux[0], frame='TOPOCENT')
    bary = Spectrum(spectral_axis=spectral_axis[0], flux=flux[0], frame='BARYCENT')
    with pytest.raises(ValueError, match="must have the same frame"), \
            warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        SpectrumCollection.from_spectra([topo, bary])
    with pytest.raises(ValueError, match="Cannot separately set frame"):
        SpectrumCollection(flux, spectral_axis=collection.spectral_axis, frame='BARYCENT')

    plain = SpectrumCollection(flux, spectral_axis=spectral_axis)
    assert plain.medium is None and plain.frame is None and plain.in_rest_frame is None
    assert plain.barycentric_correction is None
