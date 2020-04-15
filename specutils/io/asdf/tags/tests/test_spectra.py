import pytest
# Make sure these tests do not run if ASDF is not installed
pytest.importorskip('asdf')

import numpy as np

import astropy.units as u
from astropy.coordinates import FK5
from astropy.nddata import StdDevUncertainty

from asdf.tests.helpers import assert_roundtrip_tree
import asdf

from specutils import Spectrum1D, SpectrumList, SpectralAxis


def create_spectrum1d(xmin, xmax, uncertainty=None):
    flux = np.ones(xmax - xmin) * u.Jy
    wavelength = np.arange(xmin, xmax) * 0.1 * u.nm
    uncertainty = StdDevUncertainty(np.ones(xmax - xmin) * u.Jy) if uncertainty is not None else None
    return Spectrum1D(spectral_axis=wavelength, flux=flux,
                      uncertainty=uncertainty)


def test_asdf_spectrum1d(tmpdir):

    spectrum = create_spectrum1d(5100, 5300)

    tree = dict(spectrum=spectrum)
    assert_roundtrip_tree(tree, tmpdir)


def test_asdf_spectrum1d_uncertainty(tmpdir):

    spectrum = create_spectrum1d(5100, 5300, uncertainty=True)

    tree = dict(spectrum=spectrum)
    assert_roundtrip_tree(tree, tmpdir)

@pytest.mark.xfail
def test_asdf_spectralaxis(tmpdir):

    wavelengths  = np.arange(5100, 5300) * 0.1 * u.nm
    spectral_axis = SpectralAxis(wavelengths, bin_specification="edges")
    tree = dict(spectral_axis=spectral_axis)
    assert_roundtrip_tree(tree, tmpdir)

def test_asdf_spectrumlist(tmpdir):

    spectra = SpectrumList([
        create_spectrum1d(5100, 5300),
        create_spectrum1d(5000, 5500),
        create_spectrum1d(0, 100),
        create_spectrum1d(1, 5)
    ])

    tree = dict(spectra=spectra)
    assert_roundtrip_tree(tree, tmpdir)


@pytest.mark.filterwarnings("error::UserWarning")
def test_asdf_url_mapper():
    """Make sure specutils asdf extension url_mapping doesn't interfere with astropy schemas"""
    frame = FK5()
    af = asdf.AsdfFile()
    af.tree = {'frame': frame}
