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
import copy
from astropy import constants
import math

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

        self.fits_header_writers = {'linear': self._write_fits_header_linear,
                                    'matrix': self._write_fits_header_matrix,
                                    'multispec': self._write_fits_header_multispec}

    def __call__(self, pixel_indices):
        return super(Spectrum1DPolynomialWCS, self).__call__(
            pixel_indices) * self.unit

    def write_fits_header(self, header, spectral_axis=1, method='linear'):
        self.fits_header_writers[method](header, spectral_axis)

    def _write_fits_header_linear(self, header, spectral_axis=1):
        header['cdelt{0}'.format(spectral_axis)] = self.c1.value
        header['crval{0}'.format(spectral_axis)] = self.c0.value

        if self._unit is not None:
            unit_string = self.unit
            if isinstance(self.unit, u.UnitBase):
                if self.unit == u.AA:
                    unit_string = 'angstroms'
                else:
                    unit_string = self.unit.to_string()

            header['cunit{0}'.format(spectral_axis)] = unit_string

    def _write_fits_header_matrix(self, header, spectral_axis=1):
        header['cd{0}_{1}'.format(spectral_axis, spectral_axis)] = self.c1.value
        header['crval{0}'.format(spectral_axis)] = self.c0.value

        if self._unit is not None:
            unit_string = self.unit
            if isinstance(self.unit, u.UnitBase):
                if self.unit == u.AA:
                    unit_string = 'angstroms'
                else:
                    unit_string = self.unit.to_string()

            header['cunit{0}'.format(spectral_axis)] = unit_string

    # can only be implemented, when the reader is in place
    def _write_fits_header_multispec(self, header, spectral_axis=1):
        pass


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

    def get_fits_spec(self):
        func_type = 2
        order = self.degree + 1
        coefficients = [self.__getattr__('c{0}'.format(i)).value
                        for i in range(order)]
        return [func_type, order, self.pmin, self.pmax] + coefficients

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

    def get_fits_spec(self):
        func_type = 1
        order = self.degree + 1
        coefficients = [self.__getattr__('c{0}'.format(i)).value
                        for i in range(order)]
        return [func_type, order, self.pmin, self.pmax] + coefficients

class Spectrum1DIRAFBSplineWCS(BaseSpectrum1DWCS, BSplineModel):
    """
    WCS for polynomial dispersion using BSpline functions with transformation
    required for processing IRAF specification described at
    http://iraf.net/irafdocs/specwcs.php
    """

    def __init__(self, degree, npieces, y, pmin, pmax):
        from scipy.interpolate import splrep
        self.pmin = pmin
        self.pmax = pmax
        self.npieces = npieces
        self.y = y

        x = np.arange(npieces + degree)
        knots, coefficients, _ = splrep(x, y, k=degree)
        super(Spectrum1DIRAFBSplineWCS, self).__init__(degree, knots,
                                                       coefficients)

    def __call__(self, pixel_indices):
        s = (pixel_indices * 1.0 * self.npieces) / (self.pmax - self.pmin)
        return super(Spectrum1DIRAFBSplineWCS, self).__call__(s)

    def get_fits_spec(self):
        if self.degree == 1:
            func_type = 4
        elif self.degree == 3:
            func_type = 3
        else:
            raise Spectrum1DWCSFITSError("Fits spec undefined for degree = "
                                         "{0}".format(self.degree))
        return [func_type, self.npieces, self.pmin, self.pmax] + self.y


class Spectrum1DIRAFCombinationWCS(BaseSpectrum1DWCS):
    """
    WCS that combines multiple WCS using their weights, zero index and doppler
    factor. The formula used is:
    Dispersion = Sum over all WCS
        [Weight * (Zero point offset + WCS(pixels)) / (1 + doppler factor)]
    """
    def __init__(self, num_pixels, aperture=1, beam=88, aperture_low=0.0,
                 aperture_high=0.0, doppler_factor=0.0, unit=None):
        self.wcs_list = []
        self.aperture = aperture
        self.beam = beam
        self.num_pixels = num_pixels
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

    def get_fits_spec(self):
        disp_type = 2
        dispersion0 = self.__call__(np.zeros(1))[0].value
        dispersion = self.__call__(np.arange(self.num_pixels)).value
        avg_disp_delta = (dispersion[1:] - dispersion[:-1]).mean()
        spec = [self.aperture, self.beam, disp_type, dispersion0,
                avg_disp_delta, self.num_pixels, self.doppler_factor,
                self.aperture_low, self.aperture_high]
        for wcs, weight, zero_point_offset in self.wcs_list:
            spec.extend([weight, zero_point_offset])
            spec.extend(wcs.get_fits_spec())

        return " ".join(map(str, spec))


class WeightedCombinationWCS(Model):
    """
    A weighted combination WCS model. This model combines multiple WCS using a
    weight and a zero point offset. Suppose there are n WCS'es. When called with
    an input, this WCS returns the following:
            Sum over i = 1 to n
                [weight_i * (zero_point_offset_i + WCS_i(input))
    WCS can be added using the add_wcs method, along with it's weight and zero
    point offset.

    Parameters
    -----------
    wcs_list : list of callable objects, optional
        The object's wcs_list will be instantiated using wcs from this list,
        with weight as 1.0, and zero point offset as 0.0
    """
    def __init__(self, wcs_list=[]):
        self.wcs_list = []
        for wcs in wcs_list:
            self.add_WCS(wcs)

    def add_WCS(self, wcs, weight=1.0, zero_point_offset=0.0):
        """
        Add a WCS/function pointer to be evaluated when this WCS is called. The
        results of calling this WCS on the input will be added to the overall
        result, after applying the weight and the ero point offset.

        Parameters
        -----------
        wcs : callable
            The WCS to be added
        weight: float, optional
            The weight of the WCS in this model
        zero_point: float, optional
            The output of the WCS will be offset by this value before being
            multiplied by the weight
        """
        self.wcs_list.append((wcs, weight, zero_point_offset))

    def __call__(self, input):
        """
        Applies the WCS'es and functions to the input and returns the
        weighted sum of all the results

        Parameters
        -----------
        input : numpy array
            The input to the composite WCS
        """
        output = np.zeros(len(input))
        for wcs, weight, zero_point_offset in self.wcs_list:
            output += weight * (zero_point_offset + wcs(input))
        return output


class CompositeWCS(Model):
    """
    A composite WCS model. This model applies multiple WCS in-order to a
    particular input. Suppose there are 4 WCS'es, When called, this WCS returns
    the following:
                        wcs_4(wcs_3(wcs_2(wcs_1(input))))
    The WCS which is added first is applied first. WCS'es can be added using
    the add_wcs method.

    Parameters
    -----------
    wcs_list : list of callable objects, optional
        The object's wcs_list will be instantiated using wcs from this list
    """
    def __init__(self, wcs_list=[]):
        self.wcs_list = []
        for wcs in wcs_list:
            self.add_WCS(wcs)

    def add_WCS(self, wcs):
        """
        Add a WCS/function pointer on top of the current transformation. This
        will be called after the previous items in the list

        Parameters
        -----------
        wcs : callable
            The WCS to be added
        """
        self.wcs_list.append(wcs)

    def __call__(self, input):
        """
        Applies the chain of WCS'es and functions to the input and returns the
        result

        Parameters
        -----------
        input : numpy array
            The input to the composite WCS
        """
        output = input
        for wcs in self.wcs_list:
            output = wcs(output)
        return output

class DopplerShift(Model):
    """
    Applies doppler shift to the input. Returns:
                    input * doppler factor
    Parameters
    -----------
    doppler_factor : float
        the doppler factor
    """

    @classmethod
    def from_velocity(cls, velocity):
        """
        Instantiates the doppler shift model from the velocity (v) given, the
        doppler factor is computed and stored using the following formula:
                   doppler factor = sqrt((1 + v/c)/(1 - v/c))
        where c is the speed of light

        Parameters
        -----------
        velocity : float
            the relative velocity between the observer and the source
        """
        beta = velocity/constants.c
        doppler_factor = math.sqrt((1 + beta)/(1 - beta))
        return cls(doppler_factor)

    @classmethod
    def from_z(cls, z):
        """
        Instantiates the doppler shift model from redshift (z).
        """
        return cls(1 + z)

    def __init__(self, doppler_factor):
        self.doppler_factor = doppler_factor

    def __call__(self, input):
        """
        Applies the doppler shift to the input, and returns the result

        Parameters
        -----------
        input : numpy array
            The input to be shifted
        """
        return input * self.doppler_factor

    def inverse(self, input):
        """
        Applies doppler de-shift to the input, and returns the result. 

        Parameters
        -----------
        input : numpy array
            The input to be shifted
        """
        return input / self.doppler_factor

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
