import pytest
# Make sure these tests do not run if ASDF is not installed
pytest.importorskip('asdf')

import numpy as np  # noqa: E402

import astropy.units as u  # noqa: E402
from astropy.coordinates import FK5  # noqa: E402
from astropy.nddata import StdDevUncertainty  # noqa: E402

from asdf.testing.helpers import roundtrip_object  # noqa: E402
import asdf  # noqa: E402

from specutils import Spectrum1D, SpectrumList, SpectralAxis  # noqa: E402
from specutils.io.asdf.tags.spectra import Spectrum1DType, SpectrumListType  # noqa: E402


def create_spectrum1d(xmin, xmax, uncertainty=None):
    flux = np.ones(xmax - xmin) * u.Jy
    wavelength = np.arange(xmin, xmax) * 0.1 * u.nm
    uncertainty = StdDevUncertainty(np.ones(xmax - xmin) * u.Jy) if uncertainty is not None else None
    return Spectrum1D(spectral_axis=wavelength, flux=flux,
                      uncertainty=uncertainty)


@pytest.mark.filterwarnings('ignore:ASDF functionality for astropy is being moved out')
def test_asdf_spectrum1d(tmpdir):

    spectrum = create_spectrum1d(5100, 5300)
    Spectrum1DType.assert_equal(spectrum, roundtrip_object(spectrum))


@pytest.mark.filterwarnings('ignore:ASDF functionality for astropy is being moved out')
def test_asdf_spectrum1d_uncertainty(tmpdir):

    spectrum = create_spectrum1d(5100, 5300, uncertainty=True)
    Spectrum1DType.assert_equal(spectrum, roundtrip_object(spectrum))


@pytest.mark.xfail
def test_asdf_spectralaxis(tmpdir):

    wavelengths  = np.arange(5100, 5300) * 0.1 * u.nm
    spectral_axis = SpectralAxis(wavelengths, bin_specification="edges")
    # there is no implemented asdf type for SpectralAxis and no defined
    # equality comparison (assert_equal)
    assert roundtrip_object(spectral_axis) == spectral_axis


@pytest.mark.filterwarnings('ignore:ASDF functionality for astropy is being moved out')
def test_asdf_spectrumlist(tmpdir):

    spectra = SpectrumList([
        create_spectrum1d(5100, 5300),
        create_spectrum1d(5000, 5500),
        create_spectrum1d(0, 100),
        create_spectrum1d(1, 5)
    ])
    SpectrumListType.assert_equal(spectra, roundtrip_object(spectra))


@pytest.mark.filterwarnings("error::UserWarning")
def test_asdf_url_mapper():
    """Make sure specutils asdf extension url_mapping doesn't interfere with astropy schemas"""
    frame = FK5()
    af = asdf.AsdfFile()
    af.tree = {'frame': frame}
