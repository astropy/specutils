import numpy as np

import astropy.units as u

from asdf.tests.helpers import assert_roundtrip_tree

from specutils import Spectrum1D, SpectrumList


def test_asdf_spectrum1d(tmpdir):

    flux = np.random.randn(200)*u.Jy
    wavelength = np.arange(5100, 5300)*u.AA
    spec1d = Spectrum1D(spectral_axis=wavelength, flux=flux)

    tree = dict(spectrum=spec1d)
    assert_roundtrip_tree(tree, tmpdir)
