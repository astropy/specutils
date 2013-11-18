# This module provides a basic and probably temporary WCS solution until astropy has a wcs built-in

import numpy as np

from astropy.utils import misc
from astropy.modeling import Model, polynomial
from astropy.modeling.parameters import Parameter


import astropy.units as u

from astropy.utils.misc import deprecated

valid_spectral_units = [u.pix, u.km/u.s, u.m, u.Hz, u.erg]



def check_valid_unit(unit):
    if not any([unit.is_equivalent(x) for x in valid_spectral_units]):
                raise ValueError("Unit %r is not recognized as a valid spectral unit.  Valid units are: " % unit.to_string() +
                                 ", ".join([x.to_string() for x in valid_spectral_units]))

class Spectrum1DWCSError(Exception):
    pass


class Spectrum1DWCSFITSError(Spectrum1DWCSError):
    pass

class Spectrum1DWCSUnitError(Spectrum1DWCSError):
    pass

class BaseSpectrum1DWCS(Model):
    """
    Base class for a Spectrum1D WCS
    """


    n_inputs = 1
    n_outputs = 1

    _default_equivalencies = u.spectral()

    @property
    def equivalencies(self):
        """Equivalencies for spectral axes include spectral equivalencies and doppler"""
        if hasattr(self, '_equivalencies'):
            return self._equivalencies
        else:
            return self._default_equivalencies

    #@equivalencies.setter
    #def equivalencies(self, equiv):
    #    u.core._normalize_equivalencies(equiv)
    #    self._equivalencies = equiv

    @property
    def unit(self):
        if self._unit is None:
            return 1.0
        else:
            return self._unit

    @unit.setter
    def unit(self, value):
        self._unit = value


    def reset_equivalencies(self):
        """
        Reset the equivalencies to the defaults (probably u.spectral())
        """
        self._equivalencies = self._default_equivalencies

    def add_equivalency(self, new_equiv):
        """
        Add a new equivalency

        Parameters
        ----------
        new_equiv: list
            A list of equivalency mappings

        Examples
        --------
        >>> wcs.add_equivalency(u.doppler_optical(5*u.AA))
        """
        u.core._normalize_equivalencies(new_equiv)
        self._equivalencies += new_equiv


class Spectrum1DLookupWCS(BaseSpectrum1DWCS):
    """
    A simple lookup table wcs

    Parameters
    ----------

    lookup_table : ~np.ndarray or ~astropy.units.Quantity
        lookup table for the array
    """

    n_inputs = 1
    n_outputs = 1
    lookup_table_parameter = Parameter('lookup_table_parameter')

    def __init__(self, lookup_table, unit=None, lookup_table_interpolation_kind='linear'):
        super(Spectrum1DLookupWCS, self).__init__()

        if unit is not None:
            self.lookup_table_parameter = u.Quantity(lookup_table, unit)
            self.unit = lookup_table.unit
        else:
            self.lookup_table_parameter = lookup_table
            self.unit = None


        self.lookup_table_interpolation_kind = lookup_table_interpolation_kind

        #Making sure that 1d transformations are sensible
        assert self.lookup_table_parameter.value.ndim == 1

        #check that array gives a bijective transformation (that forwards and backwards transformations are unique)
        if len(self.lookup_table_parameter.value) != len(np.unique(self.lookup_table_parameter.value)):
            raise Spectrum1DWCSError('The Lookup Table does not describe a unique transformation')
        self.pixel_index = np.arange(len(self.lookup_table_parameter.value))

    def __call__(self, pixel_indices):
        if self.lookup_table_interpolation_kind == 'linear':
            return np.interp(pixel_indices, self.pixel_index, self.lookup_table_parameter.value, left=np.nan,
                             right=np.nan) * self.unit
        else:
            raise NotImplementedError('Interpolation type %s is not implemented' % self.lookup_table_interpolation_kind)


    def invert(self, dispersion_values):
        if self.lookup_table_interpolation_kind == 'linear':
            return np.interp(dispersion_values, self.lookup_table_parameter.value, self.pixel_index, left=np.nan, right=np.nan)
        else:
            raise NotImplementedError('Interpolation type %s is not implemented' % self.lookup_table_interpolation_kind)



class Spectrum1DLinearWCS(BaseSpectrum1DWCS):
    """
    A simple linear wcs
    """

    dispersion0 = Parameter('dispersion0')
    dispersion_delta = Parameter('dispersion_delta')

    @deprecated('0.dev??', message='please use Spectrum1DPolynomialWCS')
    def __init__(self, dispersion0, dispersion_delta, pixel_index, unit):
        super(Spectrum1DLinearWCS, self).__init__()

        #### Not clear what to do about units of dispersion0 and dispersion_delta.
        # dispersion0 should have units like angstrom, whereas dispersion_delta should have units like angstrom/pix
        # for now I assume pixels don't have units and both dispersion0 and dispersion_delta should have the same unit
        dispersion0 = u.Quantity(dispersion0, unit)
        dispersion_delta = u.Quantity(dispersion_delta, unit)

        check_valid_unit(dispersion0.unit)
        check_valid_unit(dispersion_delta.unit)


        ##### Quick fix - needs to be fixed in modelling ###
        if unit is None:
            unit = dispersion0.unit

        self.unit = unit

        self.dispersion0 = dispersion0.value
        self.dispersion_delta = dispersion_delta.value
        self.pixel_index = pixel_index



    def __call__(self, pixel_indices):
        if misc.isiterable(pixel_indices) and not isinstance(pixel_indices, basestring):
            pixel_indices = np.array(pixel_indices)
        return (self.dispersion0 + self.dispersion_delta * (pixel_indices - self.pixel_index)) * self.unit

    def invert(self, dispersion_values):
        if not hasattr(dispersion_values, 'unit'):
            raise u.UnitsException('Must give a dispersion value with a valid unit (i.e. quantity 5 * u.Angstrom)')

        if misc.isiterable(dispersion_values) and not isinstance(dispersion_values, basestring):
            dispersion_values = np.array(dispersion_values)
        return float((dispersion_values - self.dispersion0) / self.dispersion_delta) + self.pixel_index



class Spectrum1DPolynomialWCS(BaseSpectrum1DWCS, polynomial.Polynomial1D):

    def __init__(self, degree, unit=None, domain=None, window=[-1, 1], param_dim=1, **params):
        super(Spectrum1DPolynomialWCS, self).__init__(degree, domain=domain, window=window, param_dim=param_dim,
                 **params)
        self.unit = unit

    def __call__(self, pixel_indices):
        return polynomial.Polynomial1D.__call__(self, pixel_indices) * self.unit


class Spectrum1DLegendreWCS(BaseSpectrum1DWCS, polynomial.Legendre1D):

    def __init__(self, degree, unit=None, domain=None, window=[-1, 1], param_dim=1,
                 **params):
        super(Spectrum1DLegendreWCS, self).__init__(degree, domain=domain, window=window, param_dim=param_dim,
                 **params)
        check_valid_unit(unit)
        self.unit = unit

    def __call__(self, pixel_indices):
        return polynomial.Legendre1D.__call__(self, pixel_indices) * self.unit


# Checking which WCSes

fits_capable_wcs = []

for wcs in BaseSpectrum1DWCS.__subclasses__():
    if hasattr(wcs, 'from_fits_header') and hasattr(wcs, 'to_fits_header'):
        fits_capable_wcs.append(wcs)



def _parse_doppler_convention(dc):
    dcd = {'relativistic': u.doppler_relativistic,
           'radio': u.doppler_radio,
           'optical': u.doppler_optical}
    if dc in dcd:
        return dcd[dc]
    elif dc in dcd.values(): # allow users to specify the convention directly
        return dc
    else:
        raise ValueError("Doppler convention must be one of " + ",".join(dcd.keys()))
