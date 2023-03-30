import asdf
import numpy as np
import pytest
from astropy import units as u
from astropy.coordinates import FK5
from astropy.nddata import StdDevUncertainty

from specutils import Spectrum1D, SpectrumList, SpectralAxis
from specutils.io.asdf.tests.helpers import (
    assert_spectrum1d_equal, assert_spectrumlist_equal, assert_spectral_axis_equal)


def create_spectrum1d(xmin, xmax, uncertainty=False, mask=False):
    flux = np.ones(10) * u.Jy
    wavelength = np.linspace(xmin, xmax, 10) * u.nm
    unc = StdDevUncertainty(flux * 0.1) if uncertainty else None
    msk = np.array([0, 1, 1, 0, 1, 0, 1, 1, 0, 1], dtype=np.uint8) if mask else None
    return Spectrum1D(spectral_axis=wavelength, flux=flux, uncertainty=unc, mask=msk)


@pytest.mark.parametrize('uncertainty', [False, True])
@pytest.mark.parametrize('mask', [False, True])
def test_asdf_spectrum1d(tmp_path, uncertainty, mask):
    file_path = tmp_path / "test.asdf"
    spectrum = create_spectrum1d(510, 530, uncertainty=uncertainty, mask=mask)
    with asdf.AsdfFile() as af:
        af["spectrum"] = spectrum
        af.write_to(file_path)

    with asdf.open(file_path) as af:
        assert_spectrum1d_equal(af["spectrum"], spectrum)


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
        create_spectrum1d(510, 530),
        create_spectrum1d(500, 550),
        create_spectrum1d(0, 10),
        create_spectrum1d(0.1, 0.5)
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
