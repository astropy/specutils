import asdf
import numpy as np
import pytest
from astropy import units as u
from astropy.coordinates import FK5
from astropy.nddata import StdDevUncertainty

from specutils import Spectrum1D, SpectrumList, SpectralAxis
from specutils.io.asdf.helpers import (
    assert_spectrum1d_equal, assert_spectrumlist_equal, assert_spectral_axis_equal)


def create_spectrum1d(xmin, xmax, uncertainty=None):
    flux = np.ones(xmax - xmin) * u.Jy
    wavelength = np.arange(xmin, xmax) * u.AA
    uncertainty = StdDevUncertainty(np.ones(xmax - xmin) * u.Jy) if uncertainty is not None else None
    return Spectrum1D(spectral_axis=wavelength, flux=flux, uncertainty=uncertainty)


@pytest.mark.parametrize('uncertainty', [False, True])
@pytest.mark.filterwarnings("ignore:.*The unit 'Angstrom' has been deprecated.*")
def test_asdf_spectrum1d(tmp_path, uncertainty):
    file_path = tmp_path / "test.asdf"
    spectrum = create_spectrum1d(5100, 5300, uncertainty=uncertainty)
    with asdf.AsdfFile() as af:
        af["spectrum"] = spectrum
        af.write_to(file_path)

    with asdf.open(file_path) as af:
        assert_spectrum1d_equal(af["spectrum"], spectrum)


@pytest.mark.filterwarnings("ignore:.*The unit 'Angstrom' has been deprecated.*")
def test_asdf_spectralaxis(tmp_path):
    file_path = tmp_path / "test.asdf"
    wavelengths  = np.arange(5100, 5300) * u.AA
    spectral_axis = SpectralAxis(wavelengths, bin_specification="edges")

    with asdf.AsdfFile() as af:
        af["spectral_axis"] = spectral_axis
        af.write_to(file_path)

    with asdf.open(file_path) as af:
        assert_spectral_axis_equal(af["spectral_axis"], spectral_axis)


@pytest.mark.filterwarnings("ignore:.*The unit 'Angstrom' has been deprecated.*")
def test_asdf_spectrumlist(tmp_path):
    file_path = tmp_path / "test.asdf"
    spectra = SpectrumList([
        create_spectrum1d(5100, 5300),
        create_spectrum1d(5000, 5500),
        create_spectrum1d(0, 100),
        create_spectrum1d(1, 5)
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
