# This module provides a basic and probably temporary WCS solution until astropy has a wcs built-in
#This is all built upon @nden's work on the models class

# Add when nden's models are merged
# from astropy import models

import numpy as np

from astropy.utils import misc
from astropy.io import fits
from astropy import modeling
import astropy.units as u

valid_spectral_units = [u.pix, u.km/u.s, u.m, u.Hz, u.erg]

class BaseSpectrum1DWCSError(Exception):
    pass

class BaseSpectrum1DWCS(modeling.Model):
    """
    Base class for a Spectrum1D WCS
    """

    def pixel2dispersion(self, pixel_index):
        """
        This should be the forward transformation in normal WCS classes
        """
        return self[pixel_index]

    def dispersion2pixel(self, dispersion_value):
        """
        This should be the inverse transformation in normal WCS classes
        """
        return self.invert(dispersion_value)

    @misc.lazyproperty
    def lookup_table(self):
        return self(self.pixel_index)


class Spectrum1DLookupWCS(BaseSpectrum1DWCS):
    """
    A simple lookup table wcs

    Parameters
    ----------

    lookup_table : ~np.ndarray
        lookup table for the array
    """

    def __init__(self, lookup_table, unit=None, lookup_table_interpolation_kind='linear'):

        self._lookup_table_parameter = modeling.parameters.Parameter('lookup_table_parameter', lookup_table, self, 1)

        if unit is None and not hasattr(lookup_table, 'unit'):
            raise TypeError("Lookup table must have a unit attribute or units must be specified.")
        elif unit is not None:
            try:
                unit = u.Unit(unit)
            except ValueError:
                # e.g., Unit("blah")
                raise ValueError("Invalid unit specified")

            if not any([unit.is_equivalent(x) for x in valid_spectral_units]):
                raise ValueError("Unit %r is not recognized as a valid spectral unit.  Valid units are: " % unit.to_string() +
                                 ", ".join([x.to_string() for x in valid_spectral_units]))
            self.lookup_table = lookup_table * unit
        elif unit is None: # lookup_table is already unit'd
            self.lookup_table = lookup_table

        self.lookup_table_interpolation_kind = lookup_table_interpolation_kind
        super(Spectrum1DLookupWCS, self).__init__((), n_inputs=1, n_outputs=1, param_dim=1)

        #check that array gives a bijective transformation (that forwards and backwards transformations are unique)
        if len(self._lookup_table_parameter[0]) != len(np.unique(self._lookup_table_parameter[0])):
            raise BaseSpectrum1DWCSError('The Lookup Table does not describe a unique transformation')
        self.pixel_index = np.arange(len(self._lookup_table_parameter[0]))

    def __call__(self, pixel_indices):
        if self.lookup_table_interpolation_kind == 'linear':
            return np.interp(pixel_indices, self.pixel_index, self._lookup_table_parameter[0], left=np.nan, right=np.nan)
        else:
            raise NotImplementedError('Interpolation type %s is not implemented' % self.lookup_table_interpolation_kind)


    def invert(self, dispersion_values):
        if self.lookup_table_interpolation_kind == 'linear':
            return np.interp(dispersion_values, self._lookup_table_parameter[0], self.pixel_index, left=np.nan, right=np.nan)
        else:
            raise NotImplementedError('Interpolation type %s is not implemented' % self.lookup_table_interpolation_kind)



class Spectrum1DLinearWCS(BaseSpectrum1DWCS):
    """
    A simple linear wcs
    """

    param_names = ['dispersion0', 'dispersion_delta']

    @classmethod
    def from_header(cls, header, unit=None, spectroscopic_axis_number=1,
            **kwargs):
        """
        Load a simple linear Spectral WCS from a FITS header.

        Respects the following keywords:

         * CDELT or CD: delta-dispersion
         * CRVAL: reference position
         * CRPIX: reference pixel
         * CUNIT: dispersion units (parsed by astropy.units)

        Parameters
        ----------
        header : str
            FITS Filename or FITS Header instance
        unit : astropy.units.Unit
            Specified units. Overrides CUNIT1 if specified.

        .. todo:: Allow FITS files or just headers to be passed...
        """
        if isinstance(header, basestring):
            header = fits.getheader(header, **kwargs)

        if unit is None:
            try:
                unit = u.Unit(header.get('CUNIT%i' % spectroscopic_axis_number))
            except u.UnitsException:
                raise u.UnitsException("No units were specified and CUNIT did not contain unit information.")

        cdelt = header.get('CDELT%i' % spectroscopic_axis_number)
        if cdelt is None:
            cdelt = header.get('CD%i_%i' % (spectroscopic_axis_number,spectroscopic_axis_number))

        crpix = header['CRPIX%i' % spectroscopic_axis_number]
        crval = header['CRVAL%i' % spectroscopic_axis_number]

        return cls(crval, cdelt, crpix - 1, unit=unit)


    def __init__(self, dispersion0, dispersion_delta, pixel_index, unit=None):
        self.unit = unit
        self.dispersion0 = dispersion0 * self.unit
        self.dispersion_delta = dispersion_delta * self.unit
        self.pixel_index = pixel_index

    def __call__(self, pixel_indices):
        if misc.isiterable(pixel_indices) and not isinstance(pixel_indices, basestring):
            pixel_indices = np.array(pixel_indices)
        return self.dispersion0 + self.dispersion_delta * (pixel_indices - self.pixel_index)

    def invert(self, dispersion_values):
        if not hasattr(dispersion_values, 'unit'):
            raise u.UnitsException('Must give a dispersion value with a valid unit')
        if misc.isiterable(dispersion_values) and not isinstance(dispersion_values, basestring):
            dispersion_values = np.array(dispersion_values)
        return float((dispersion_values - self.dispersion0) / self.dispersion_delta) + self.pixel_index


#### EXAMPLE implementation for Chebyshev
#class ChebyshevSpectrum1D(models.ChebyshevModel):
class ChebyshevSpectrum1D(BaseSpectrum1DWCS):

    @classmethod
    def from_fits_header(cls, header):
        pass
        ### here be @hamogu's code ###
        #degree, parameters = hamogu_read_fits(header)
        #return cls(degree, **parameters)







