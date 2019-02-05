import numpy as np

import astropy.units as u

from asdf.tests.helpers import assert_roundtrip_tree

from specutils import Spectrum1D, SpectrumList


def create_spectrum1d(xmin, xmax, uncertainty=None):
    flux = np.ones(xmax-xmin) * u.Jy
    wavelength = np.arange(xmin, xmax) * u.AA
    uncertainty = np.ones(xmax-xmin) if uncertainty is not None else None
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


def test_asdf_spectrumlist(tmpdir):

    spectra = SpectrumList([
        create_spectrum1d(5100, 5300),
        create_spectrum1d(5000, 5500),
        create_spectrum1d(0, 100),
        create_spectrum1d(1, 5)
    ])

    tree = dict(spectra=spectra)
    assert_roundtrip_tree(tree, tmpdir)
