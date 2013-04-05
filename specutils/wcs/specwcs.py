# This module provides a basic and probably temporary WCS solution until astropy has a wcs built-in
#This is all built upon @nden's work on the models class

# Add when nden's models are merged
# from astropy import models

import numpy as np

from astropy.utils import misc
from astropy.io import fits

class NDenModelsPlaceHolder(object):

    pass
class BaseSpectrum1DWCSError(Exception):
    pass

class BaseSpectrum1DWCS(NDenModelsPlaceHolder):
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

    def create_lookup_table(self, pixel_indices):
        self.lookup_table = self(pixel_indices)




class Spectrum1DLookupWCS(BaseSpectrum1DWCS):
    """
    A simple lookup table wcs

    Parameters
    ----------

    lookup_table : ~np.ndarray
        lookup table for the array
    """

    def __init__(self, lookup_table, unit=None):
        self.unit = unit
        self.lookup_table = lookup_table

        #check that array gives a bijective transformation (that forwards and backwards transformations are unique)

        if len(self.lookup_table) != len(np.unique(self.lookup_table)):
            raise BaseSpectrum1DWCSError('The Lookup Table does not describe a unique transformation')


    def __call__(self, pixel_indices):
        if misc.isiterable(pixel_indices) and not isinstance(pixel_indices, basestring):
            pixel_indices = np.array(pixel_indices)
        return self.lookup_table[pixel_indices]

    def invert(self, dispersion_values):

        return np.searchsorted(self.lookup_table, dispersion_values)



class Spectrum1DLinearWCS(BaseSpectrum1DWCS):
    """
        A simple linear wcs

    """


    @classmethod
    def from_fits(cls, fname, unit=None, **kwargs):
        header = fits.getheader(fname, **kwargs)
        return cls(header['CRVAL1'], header['CDELT1'], header['CRPIX1'] - 1, unit=unit)


    def __init__(self, dispersion0, dispersion_delta, dispersion_pixel0=0, unit=None):
        self.unit = unit
        self.dispersion0 = dispersion0
        self.dispersion_delta = dispersion_delta
        self.dispersion_pixel0 = dispersion_pixel0

    def __call__(self, pixel_indices):
        if misc.isiterable(pixel_indices) and not isinstance(pixel_indices, basestring):
            pixel_indices = np.array(pixel_indices)
        return self.dispersion0 + self.dispersion_delta * (pixel_indices - self.dispersion_pixel0)

    def invert(self, dispersion_values):
        if misc.isiterable(dispersion_values) and not isinstance(dispersion_values, basestring):
            dispersion_values = np.array(dispersion_values)
        return (dispersion_values - self.dispersion0) / self.dispersion_delta + self.dispersion_pixel0


#### EXAMPLE implementation for Chebyshev
#class ChebyshevSpectrum1D(models.ChebyshevModel):
class ChebyshevSpectrum1D(BaseSpectrum1DWCS):

    @classmethod
    def from_fits_header(cls, header):
        pass
        ### here be @hamogu's code ###
        #degree, parameters = hamogu_read_fits(header)
        #return cls(degree, **parameters)







