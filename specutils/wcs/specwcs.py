# This module provides a basic and probably temporary WCS solution until astropy has a wcs built-in
#This is all built upon @nden's work on the models class
from astropy import models
import numpy as np

class BaseSpectrum1DWCSError(Exception):
    pass

class BaseSpectrum1DWCS(models.Model):
    """
        Base class for a Spectrum1D WCS
    """

    def pixel2dispersion(self, pixel_index):
        """
            This should be the forward transformation in normal WCS classes
        """
        return self(pixel_index)

    def dispersion2pixel(self, dispersion_value):
        """
            This should be the inverse transformation in normal WCS classes
        """
        return self.invert(dispersion_value)



class LookupWCS(BaseSpectrum1DWCS):
    """
    A simple lookup table wcs

    Parameters
    ----------

    lookup_table : ~np.ndarray
        lookup table for
    """

    def __init__(self, lookup_table):
        self.lookup_table = lookup_table

        #check that array gives a bijective transformation (that forwards and backwards transformations are unique)

        if len(self.lookup) != len(np.unique(self.lookup_table)):
            raise BaseSpectrum1DWCSError('The Lookup Table does not describe a unique transformation')


    def __call__(self, pixel_index):
        return self.lookup_table[pixel_index]

    def invert(self, dispersion_values):
        return np.searchsorted(self.lookup_table, dispersion_values)

    @property
    def lut(self):
        return self.lookup_table


class ChebyshevSpectrum1D(models.ChebyshevModel):

    @classmethod
    def from_fits_header(cls, header):
        ### here be @hamogu's code ###
        degree, parameters = hamogu_read_fits(header)
        return cls(degree, **parameters)

    def create_lookup_table(self, pixel_indices):
        self.lut = self(pixel_indices)






