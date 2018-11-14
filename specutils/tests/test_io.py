# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module tests SpecUtils io routines
"""

from ..io.parsing_utils import generic_spectrum_from_table # or something like that
from astropy.table import Table
from astropy.utils.exceptions import AstropyUserWarning
from astropy.tests.helper import catch_warnings
import astropy.units as u
import numpy as np
import pytest
import warnings

def test_generic_spectrum_from_table(recwarn):
   """ 
   Read a simple table with wavelength, flux and uncertainty
   """
   # Create a small data set, first without uncertainties
   wave = np.arange(1,1.1,0.01)*u.AA
   flux = np.ones(len(wave))*1.e-14*u.Jy
   table = Table([wave,flux],names=["wave","flux"])

   # Test that the units and values of the Spectrum1D object match those in the table
   spectrum = generic_spectrum_from_table(table)
   assert spectrum.spectral_axis.unit == table['wave'].unit
   assert spectrum.flux.unit == table['flux'].unit
   assert spectrum.spectral_axis.unit == table['wave'].unit
   assert np.alltrue(spectrum.spectral_axis == table['wave'])
   assert np.alltrue(spectrum.flux == table['flux'])

   # Add uncertainties and retest 
   err = 0.01*flux
   table = Table([wave,flux,err],names=["wave","flux","err"])
   spectrum = generic_spectrum_from_table(table)
   assert spectrum.spectral_axis.unit == table['wave'].unit
   assert spectrum.flux.unit == table['flux'].unit
   assert spectrum.uncertainty.unit == table['err'].unit
   assert spectrum.spectral_axis.unit == table['wave'].unit
   assert np.alltrue(spectrum.spectral_axis == table['wave'])
   assert np.alltrue(spectrum.flux == table['flux'])
   assert np.alltrue(spectrum.uncertainty.array == table['err'])

   # Test for warning if standard deviation is zero or negative
   err[0] = 0.
   table = Table([wave,flux,err],names=["wave","flux","err"])
   spectrum = generic_spectrum_from_table(table)
   assert len(recwarn) == 1
   w = recwarn.pop(AstropyUserWarning)
   assert "Standard Deviation has values of 0 or less" in str(w.message) 

   # Test that exceptions are raised if there are no units
   flux = np.ones(len(wave))*1.e-14
   table = Table([wave,flux],names=["wave","flux"])
   with pytest.raises(IOError) as exc:
       spectrum = generic_spectrum_from_table(table)
       assert 'Could not identify column containing the flux' in exc
   wave = np.arange(1,1.1,0.01)
   table = Table([wave,flux,err],names=["wave","flux","err"])
   with pytest.raises(IOError) as exc:
       spectrum = generic_spectrum_from_table(table)
       assert 'Could not identify column containing the wavelength, frequency or energy' in exc
