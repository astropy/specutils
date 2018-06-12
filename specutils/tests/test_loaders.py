from __future__ import absolute_import, division

import pytest
from astropy.utils.data import get_pkg_data_filename

from ..spectra.spectrum1d import Spectrum1D


def test_spectrum1d_GMOSfits():
    optical_fits_file = get_pkg_data_filename('data/L5g_0355+11_Cruz09.fits')
    optical_spec_2 = Spectrum1D.read(optical_fits_file, format='wcs1d-fits')

    assert len(optical_spec_2.data) == 3020


def test_specific_spec_axis_unit():
    optical_fits_file = get_pkg_data_filename('data/L5g_0355+11_Cruz09.fits')
    optical_spec = Spectrum1D.read(optical_fits_file,
                                   spectral_axis_unit="Angstrom",
                                   format='wcs1d-fits')

    assert optical_spec.spectral_axis.unit == "Angstrom"
