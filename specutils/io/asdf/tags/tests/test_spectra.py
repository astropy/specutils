import asdf
import numpy as np
import pytest
from asdf.testing.helpers import roundtrip_object
from astropy import units as u
from astropy.coordinates import FK5
from astropy.nddata import StdDevUncertainty

from specutils import Spectrum1D, SpectrumList, SpectralAxis
from specutils.io.asdf.tags.tests.helpers import (
    assert_spectrum1d_equal, assert_spectrumlist_equal)


def create_spectrum1d(xmin, xmax, uncertainty=None):
    flux = np.ones(xmax - xmin) * u.Jy
    wavelength = np.arange(xmin, xmax) * u.AA
    uncertainty = StdDevUncertainty(np.ones(xmax - xmin) * u.Jy) if uncertainty is not None else None
    return Spectrum1D(spectral_axis=wavelength, flux=flux, uncertainty=uncertainty)


@pytest.mark.parametrize('uncertainty', [False, True])
def test_asdf_spectrum1d(uncertainty):
    spectrum = create_spectrum1d(5100, 5300, uncertainty=uncertainty)
    assert_spectrum1d_equal(spectrum, roundtrip_object(spectrum))


@pytest.mark.xfail
def test_asdf_spectralaxis():
    wavelengths  = np.arange(5100, 5300) * u.AA
    spectral_axis = SpectralAxis(wavelengths, bin_specification="edges")
    # There is no implemented asdf type for SpectralAxis and no defined
    # equality comparison (assert_spectralaxis_equal)
    # per the comment https://github.com/astropy/specutils/pull/645#issuecomment-614271632 ;
    # the issue is that SpectralAxis roundtrips as SpectralCoord
    # so the types differ.
    rt = roundtrip_object(spectral_axis)
    assert type(rt) == type(spectral_axis)


def test_asdf_spectrumlist():
    spectra = SpectrumList([
        create_spectrum1d(5100, 5300),
        create_spectrum1d(5000, 5500),
        create_spectrum1d(0, 100),
        create_spectrum1d(1, 5)
    ])
    assert_spectrumlist_equal(spectra, roundtrip_object(spectra))


def test_asdf_url_mapper():
    """Make sure specutils ASDF extension url_mapping does not interfere with astropy schemas."""
    af = asdf.AsdfFile()
    af.tree = {'frame': FK5()}
