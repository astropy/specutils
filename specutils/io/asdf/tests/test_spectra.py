import asdf
import numpy as np
import pytest
from astropy import units as u
from astropy.coordinates import FK5
from astropy.nddata import StdDevUncertainty
from astropy.tests.helper import assert_quantity_allclose

from specutils import Spectrum, SpectrumList, SpectralAxis
from specutils.io.asdf.tests.helpers import (
    assert_spectrum_equal, assert_spectrumlist_equal, assert_spectral_axis_equal)


def create_spectrum(xmin, xmax, uncertainty=False, mask=False):
    flux = np.ones(10) * u.Jy
    wavelength = np.linspace(xmin, xmax, 10) * u.nm
    unc = StdDevUncertainty(flux * 0.1) if uncertainty else None
    msk = np.array([0, 1, 1, 0, 1, 0, 1, 1, 0, 1], dtype=np.uint8) if mask else None
    return Spectrum(spectral_axis=wavelength, flux=flux, uncertainty=unc, mask=msk)


@pytest.mark.parametrize('uncertainty', [False, True])
@pytest.mark.parametrize('mask', [False, True])
def test_asdf_spectrum(tmp_path, uncertainty, mask):
    file_path = tmp_path / "test.asdf"
    spectrum = create_spectrum(510, 530, uncertainty=uncertainty, mask=mask)
    with asdf.AsdfFile() as af:
        af["spectrum"] = spectrum
        af.write_to(file_path)

    with asdf.open(file_path) as af:
        assert_spectrum_equal(af["spectrum"], spectrum)


def test_asdf_spectralaxis(tmp_path):
    file_path = tmp_path / "test.asdf"
    wavelengths  = np.arange(510, 530) * u.nm
    spectral_axis = SpectralAxis(wavelengths, bin_specification="edges")

    with asdf.AsdfFile() as af:
        af["spectral_axis"] = spectral_axis
        af.write_to(file_path)

    with asdf.open(file_path) as af:
        assert_spectral_axis_equal(af["spectral_axis"], spectral_axis)


def test_asdf_spectrumlist(tmp_path):
    file_path = tmp_path / "test.asdf"
    spectra = SpectrumList([
        create_spectrum(510, 530),
        create_spectrum(500, 550),
        create_spectrum(0, 10),
        create_spectrum(0.1, 0.5)
    ])
    with asdf.AsdfFile() as af:
        af["spectrum_list"] = spectra
        af.write_to(file_path)

    with asdf.open(file_path) as af:
        assert_spectrumlist_equal(af["spectrum_list"], spectra)


def test_asdf_url_mapper():
    """Make sure specutils ASDF extension url_mapping does not interfere with astropy schemas."""
    with asdf.AsdfFile() as af:
        af.tree = {'frame': FK5()}


def test_asdf_spectral_axis_metadata(tmp_path):
    from astropy.coordinates import EarthLocation, SkyCoord
    from astropy.time import Time
    from specutils.spectra.spectral_frame import SpectralMedium

    apo = EarthLocation(lat=32.78 * u.deg, lon=-105.82 * u.deg, height=2788 * u.m)
    obstime = Time('2024-03-01T05:00:00')
    target = SkyCoord(ra=120 * u.deg, dec=-30 * u.deg, radial_velocity=50 * u.km / u.s,
                      distance=100 * u.pc)
    medium = SpectralMedium('air', refraction_method='Ciddor1996', co2=400,
                            temperature=10 * u.deg_C, pressure=700 * u.hPa, humidity=0.2)
    wavelength = np.linspace(510, 530, 10) * u.nm

    # Observer and target: the radial velocity is derived and the observer rebuilt
    derived = SpectralAxis(wavelength, medium=medium, frame='TOPOCENT', obstime=obstime,
                           location=apo, target=target, doppler_rest=520 * u.nm,
                           doppler_convention='optical')
    # Manual radial velocity, nothing else
    manual = SpectralAxis(wavelength, frame='BARYCENT', radial_velocity=12 * u.km / u.s)
    # Observer that cannot be rebuilt (no location) is stored explicitly
    observer = SpectralAxis(wavelength, frame='TOPOCENT', observer=apo.get_gcrs(obstime),
                            target=target)

    file_path = tmp_path / "test.asdf"
    with asdf.AsdfFile() as af:
        af["derived"] = derived
        af["manual"] = manual
        af["observer"] = observer
        af["spectrum"] = Spectrum(spectral_axis=derived, flux=np.ones(10) * u.Jy)
        af.write_to(file_path)

    with asdf.open(file_path) as af:
        assert_spectral_axis_equal(af["derived"], derived)
        assert af["derived"].observer is not None
        assert_spectral_axis_equal(af["manual"], manual)
        assert_spectral_axis_equal(af["observer"], observer)
        assert_spectrum_equal(af["spectrum"], af["spectrum"])
        assert af["spectrum"].medium == medium
        assert af["spectrum"].frame == 'TOPOCENT'
        assert_quantity_allclose(af["spectrum"].radial_velocity, derived.radial_velocity,
                                 atol=1e-6 * u.km / u.s)
