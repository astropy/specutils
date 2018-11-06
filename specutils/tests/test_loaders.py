from __future__ import absolute_import, division

import numpy as np
from astropy.table import Table
import astropy.units as u
import pytest
from astropy.utils.data import get_pkg_data_filename

from .conftest import get_remote_data
from .. import Spectrum1D


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

def test_generic_ecsv_reader(tmpdir):
   # Create a small data set
   wave = np.arange(1,1.1,0.01)*u.AA
   flux = np.ones(len(wave))*1.e-14*u.Jy
   uncertainty = 0.01*flux
   table = Table([wave,flux,uncertainty],names=["wave","flux","uncertainty"])
   tmpfile = str(tmpdir.join('_tst.ecsv'))
   table.write(tmpfile,format='ascii.ecsv')

   # Read it in and check against the original
   spectrum = Spectrum1D.read(tmpfile,format='generic-ecsv')
   assert spectrum.spectral_axis.unit == table['wave'].unit
   assert spectrum.flux.unit == table['flux'].unit
   assert spectrum.uncertainty.unit == table['uncertainty'].unit
   assert spectrum.spectral_axis.unit == table['wave'].unit
   assert np.alltrue(spectrum.spectral_axis == table['wave'])
   assert np.alltrue(spectrum.flux == table['flux'])
   assert np.alltrue(spectrum.uncertainty.array == table['uncertainty'])



