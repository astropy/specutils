import logging

import numpy as np
from astropy import units as u
from .spectrum1d import Spectrum1D

__all__ = ['XraySpectrum1D']

# For dealing with varied unit string choices
EV   = ['eV', 'ev']
KEV  = ['kev', 'keV']
ANGS = ['angs', 'Angs', 'Angstrom', 'angstrom', 'Angstroms', 'angstroms', 'A', 'a']

def _unit_parser(unit_string):
    if unit_string in EV:
        return u.eV
    if unit_string in KEV:
        return u.keV
    if unit_string in ANGS:
        return u.angstrom

class XraySpectrum1D(Spectrum1D):
    """
    Spectrum container for holding X-ray spectrum
    WIP by eblur

    Parameters
    ----------
    bin_lo
    bin_hi
    bin_unit
    counts
    arf
    rmf
    """
    def __init__(self, bin_lo, bin_hi, bin_unit, counts, arf=None, rmf=None):
        try:
            axis_unit = u.Unit(bin_unit)
        except:
            axis_unit = _unit_parser(bin_unit)

        bin_mid = 0.5 * (bin_lo + bin_hi) * axis_unit
        Spectrum1D.__init__(self, spectral_axis=bin_mid, flux=counts)

        self.arf = arf
        self.rmf = rmf
        return

    # Convenience function for Xray people
    @property
    def counts(self):
        return self.flux
