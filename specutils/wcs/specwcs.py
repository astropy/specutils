# This module provides a basic and probably temporary WCS solution until astropy has a wcs built-in

import warnings
import numpy as np

from astropy.extern import six
from astropy.utils import misc
from astropy.modeling import Model, polynomial
from astropy.modeling.parameters import Parameter
from ..models.BSplineModel import BSplineModel
import astropy.units as u

from astropy.utils.misc import deprecated
from astropy.utils import OrderedDict
from astropy.io import fits


##### Delete at earliest convenience (currently deprecated)
#### VVVVVVVVVV
valid_spectral_units = [u.pix, u.km / u.s, u.m, u.Hz, u.erg]

@deprecated('0.dev???', 'using no units is now allowed for WCS')
def check_valid_unit(unit):
    if not any([unit.is_equivalent(x) for x in valid_spectral_units]):
        raise ValueError(
            "Unit %r is not recognized as a valid spectral unit.  Valid units are: " % unit.to_string() +
            ", ".join([x.to_string() for x in valid_spectral_units]))


#^^^^^^^^^^^^^^^^^
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
        if value is None:
            warnings.warn(
                'Initializing a Spectrum1D WCS with units set to `None` is not recommended')
            self._unit = None
        else:
            self._unit = u.Unit(value)


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

    def __init__(self, lookup_table, unit=None,
                 lookup_table_interpolation_kind='linear'):
        super(Spectrum1DLookupWCS, self).__init__()

        if unit is not None:
            self.lookup_table_parameter = u.Quantity(lookup_table, unit)
            self.unit = u.Unit(unit)
        else:
            self.lookup_table_parameter = lookup_table
            self.unit = None

        self.lookup_table_interpolation_kind = lookup_table_interpolation_kind

        #Making sure that 1d transformations are sensible
        assert self.lookup_table_parameter.value.ndim == 1

        #check that array gives a bijective transformation (that forwards and backwards transformations are unique)
        if len(self.lookup_table_parameter.value) != len(
                np.unique(self.lookup_table_parameter.value)):
            raise Spectrum1DWCSError(
                'The Lookup Table does not describe a unique transformation')
        self.pixel_index = np.arange(len(self.lookup_table_parameter.value))

    def __call__(self, pixel_indices):
        if self.lookup_table_interpolation_kind == 'linear':
            return np.interp(pixel_indices, self.pixel_index,
                             self.lookup_table_parameter.value, left=np.nan,
                             right=np.nan) * self.unit
        else:
            raise NotImplementedError(
                'Interpolation type %s is not implemented' % self.lookup_table_interpolation_kind)


    def invert(self, dispersion_values):
        if self.lookup_table_interpolation_kind == 'linear':
            return np.interp(dispersion_values,
                             self.lookup_table_parameter.value,
                             self.pixel_index, left=np.nan,
                             right=np.nan)
        else:
            raise NotImplementedError(
                'Interpolation type %s is not implemented' % self.lookup_table_interpolation_kind)


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
        if misc.isiterable(pixel_indices) and not isinstance(pixel_indices,
                                                             six.string_types):
            pixel_indices = np.array(pixel_indices)
        return (self.dispersion0 + self.dispersion_delta * (
            pixel_indices - self.pixel_index)) * self.unit

    def invert(self, dispersion_values):
        if not hasattr(dispersion_values, 'unit'):
            raise u.UnitsException(
                'Must give a dispersion value with a valid unit (i.e. quantity 5 * u.Angstrom)')

        if misc.isiterable(dispersion_values) and not isinstance(
                dispersion_values, six.string_types):
            dispersion_values = np.array(dispersion_values)
        return float((
                         dispersion_values - self.dispersion0) / self.dispersion_delta) + self.pixel_index


class Spectrum1DPolynomialWCS(BaseSpectrum1DWCS, polynomial.Polynomial1D):
    """
    WCS for polynomial dispersion. The only added parameter is a unit,
    otherwise the same as 'astropy.modeling.polynomial.Polynomial1D'
    """
    def __init__(self, degree, unit=None, domain=None, window=[-1, 1],
                 **params):
        super(Spectrum1DPolynomialWCS, self).__init__(degree, domain=domain,
                                                      window=window, **params)
        self.unit = unit

    def __call__(self, pixel_indices):
        return super(Spectrum1DPolynomialWCS, self).__call__(
            pixel_indices) * self.unit

    def add_to_header(self, header, spectral_axis=1):
        header['cdelt{0}'.format(spectral_axis)] = self.c1.value
        header['cd{0}_{1}'.format(spectral_axis, spectral_axis)] = self.c1.value
        header['crval{0}'.format(spectral_axis)] = self.c0.value
        if self._unit is not None:
            unit_string = self.unit
            if isinstance(self.unit, u.UnitBase):
                if self.unit == u.AA:
                    unit_string = 'angstroms'
                else:
                    unit_string = self.unit.name

            header['crunit{0}'.format(spectral_axis)] = unit_string


class Spectrum1DIRAFLegendreWCS(BaseSpectrum1DWCS, polynomial.Legendre1D):
    """
    WCS for polynomial dispersion using Legendre polynomials with
    transformation required for processing IRAF specification described at
    http://iraf.net/irafdocs/specwcs.php
    """
    def __init__(self, order, pmin, pmax, **coefficients):
        super(Spectrum1DIRAFLegendreWCS, self).__init__(
            order-1, domain=[pmin, pmax], **coefficients)
        self.pmin = pmin
        self.pmax = pmax

    def __call__(self, pixel_indices):
        transformed = pixel_indices + self.pmin
        return super(Spectrum1DIRAFLegendreWCS, self).__call__(transformed)


class Spectrum1DIRAFChebyshevWCS(BaseSpectrum1DWCS, polynomial.Chebyshev1D):
    """
    WCS for polynomial dispersion using Chebyshev polynomials with
    transformation required for processing IRAF specification described at
    http://iraf.net/irafdocs/specwcs.php

    See Also
    --------
    astropy.modeling.polynomial.Chebyshev1D
    astropy.modeling.polynomial.Polynomial1D
    """

    def __init__(self, order, pmin, pmax, **coefficients):
        super(Spectrum1DIRAFChebyshevWCS, self).__init__(
            order-1, domain=[pmin, pmax], **coefficients)
        self.pmin = pmin
        self.pmax = pmax

    def __call__(self, pixel_indices):
        transformed = pixel_indices + self.pmin
        return super(Spectrum1DIRAFChebyshevWCS, self).__call__(transformed)


class Spectrum1DIRAFBSplineWCS(BaseSpectrum1DWCS, BSplineModel):
    """
    WCS for polynomial dispersion using BSpline functions with transformation
    required for processing IRAF specification described at
    http://iraf.net/irafdocs/specwcs.php
    """

    @classmethod
    def from_data(cls, degree, x, y, pmin, pmax):
        from scipy.interpolate import splrep
        knots, coefficients, _ = splrep(x, y, k=degree)
        return cls(degree, knots, coefficients, pmin, pmax)

    def __init__(self, degree, knots, coefficients, pmin, pmax):
        super(Spectrum1DIRAFBSplineWCS, self).__init__(degree, knots, coefficients)
        self.pmin = pmin
        self.pmax = pmax


    def __call__(self, pixel_indices):
        n_pieces = self.n_pieces - self.degree - 2
        s = (pixel_indices * 1.0 * n_pieces) / (self.pmax - self.pmin)
        return super(Spectrum1DIRAFBSplineWCS, self).__call__(s)


class Spectrum1DIRAFCombinationWCS(BaseSpectrum1DWCS):
    """
    WCS that combines multiple WCS using their weights, zero index and doppler
    factor. The formula used is:
    Dispersion = Sum over all WCS
        [Weight * (Zero point offset + WCS(pixels)) / (1 + doppler factor)]
    """
    def __init__(self, aperture=1, beam=88, aperture_low=0.0, aperture_high=0.0,
                 doppler_factor=0.0, unit=None):
        self.wcs_list = []
        self.aperture = aperture
        self.beam = beam
        self.aperture_low = aperture_low
        self.aperture_high = aperture_high
        self.doppler_factor = doppler_factor
        self.unit = unit

    def add_WCS(self, wcs, weight=1.0, zero_point_offset=0.0):
        self.wcs_list.append((wcs, weight, zero_point_offset))

    def __call__(self, pixel_indices):
        final_dispersion = np.zeros(len(pixel_indices))
        for wcs, weight, zero_point_offset in self.wcs_list:
            dispersion = weight * (zero_point_offset + wcs(pixel_indices))
            final_dispersion += dispersion / (1 + self.doppler_factor)
        return final_dispersion * self.unit
    # Computing dispersion0 and avg dispersion delta: (for writing)
    # x2 = specx.wcs(dic['pmin'])
    # all = specx.wcs(np.arange(dic['pmin'], dic['pmax']+1))
    # y2 = (all[1:] - all[:-1]).mean()

@deprecated('0.dev???')
def _parse_doppler_convention(dc):
    dcd = {'relativistic': u.doppler_relativistic,
           'radio': u.doppler_radio,
           'optical': u.doppler_optical}
    if dc in dcd:
        return dcd[dc]
    elif dc in dcd.values():  # allow users to specify the convention directly
        return dc
    else:
        raise ValueError(
            "Doppler convention must be one of " + ",".join(dcd.keys()))
