# This module provides a basic and probably temporary WCS solution until astropy has a wcs built-in

import warnings
import numpy as np

from astropy.extern import six
from astropy.utils import misc
from astropy.modeling import Model, polynomial
from astropy.modeling.parameters import Parameter
from ..models.BSplineModel import BSplineModel
import astropy.units as u

from astropy.utils import OrderedDict
from astropy.io import fits
import copy
from astropy import constants
import math
import itertools

##### Delete at earliest convenience (currently deprecated)
#### VVVVVVVVVV
valid_spectral_units = [u.pix, u.km / u.s, u.m, u.Hz, u.erg]


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

    def __init__(self, *args, **kwargs):
        super(BaseSpectrum1DWCS, self).__init__(*args, **kwargs)

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
        super(Spectrum1DLookupWCS, self).__init__(
            lookup_table_parameter=lookup_table)

        if unit is not None:
            self.lookup_table_parameter = u.Quantity(lookup_table, unit)
            self.unit = u.Unit(unit)
        else:
            self.lookup_table_parameter = lookup_table
            self.unit = None

        self.lookup_table_interpolation_kind = lookup_table_interpolation_kind

        #Making sure that 1d transformations are sensible
        assert self.lookup_table_parameter.value.ndim == 1

        #check that array gives a bijective transformation (that forwards and
        # backwards transformations are unique)
        if len(self.lookup_table_parameter.value) != len(
                np.unique(self.lookup_table_parameter.value)):
            raise Spectrum1DWCSError(
                'The lookup table does not describe a unique transformation')
        self.pixel_index = np.arange(len(self.lookup_table_parameter.value))

    def evaluate(self, pixel_indices, lookup_table):
        lookup_table = np.squeeze(lookup_table)
        pixel_index = np.arange(len(lookup_table))
        if self.lookup_table_interpolation_kind == 'linear':
            return np.interp(pixel_indices, pixel_index, lookup_table,
                             left=np.nan, right=np.nan) * self.unit
        else:
            raise NotImplementedError(
                'Interpolation type {0} is not implemented'.format(
                    self.lookup_table_interpolation_kindkind))

    def invert(self, dispersion_values):
        if self.lookup_table_interpolation_kind == 'linear':
            return np.interp(dispersion_values,
                             self.lookup_table_parameter.value,
                             self.pixel_index, left=np.nan,
                             right=np.nan)
        else:
            raise NotImplementedError(
                'Interpolation type %s is not implemented' % self.lookup_table_interpolation_kind)


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
                                    'matrix': self._write_fits_header_matrix}


    def evaluate(self, pixel_indices, *coeffs):
        return super(Spectrum1DPolynomialWCS, self).evaluate(pixel_indices,
                                                            *coeffs) * self.unit

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
        # if abs(self.indexer.step) != 1:
        #     raise Spectrum1DWCSFITSError("WCS with indexer step greater than 1"
        #                                  "cannot be written")
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
        # if abs(self.indexer.step) != 1:
        #     raise Spectrum1DWCSFITSError("WCS with indexer step greater than 1"
        #                                  "cannot be written")
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
        self.npieces = npieces
        self.y = y

        x = np.arange(npieces + degree)
        knots, coefficients, _ = splrep(x, y, k=degree)
        super(Spectrum1DIRAFBSplineWCS, self).__init__(degree, knots,
                                                       coefficients)
        self.pmin = pmin
        self.pmax = pmax

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

    def __init__(self, wcs_list=None):
        self.wcs_list = []
        if wcs_list is not None:
            for wcs in wcs_list:
                self.add_WCS(wcs)

    def add_WCS(self, wcs, weight=1.0, zero_point_offset=0.0):
        """
        Add a WCS/function pointer to be evaluated when this WCS is called. The
        results of calling this WCS on the input will be added to the overall
        result, after applying the weight and the zero point offset.

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
        return self.__class__.evaluate(input, self.wcs_list)

    @classmethod
    def evaluate(cls, input, wcs_list):
        output = np.zeros(len(input) if hasattr(input, "__len__") else 1)
        for wcs, weight, zero_point_offset in wcs_list:
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
    def __init__(self, wcs_list=None):
        self.wcs_list = []
        if wcs_list is not None:
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
        return self.__class__.evaluate(input, self.wcs_list)

    @classmethod
    def evaluate(cls, input, wcs_list):
        output = input
        for wcs in wcs_list:
            output = wcs(output)
        return output


class DopplerShift(Model):
    """
    This model applies doppler shift to the input. Instantiates the doppler
    shift model from the velocity (v). The doppler factor is computed
    using the following formula:
                    doppler factor = sqrt((1 + v/c)/(1 - v/c))
    where c is the speed of light
    When the model is called on an input, input * doppler_factor is returned.
    The inverse of this model can also be called, which divides the input by
    the doppler factor.

    Parameters
    -----------
    velocity : float
        the relative velocity between the observer and the source
    """


    @classmethod
    def from_redshift(cls, z):
        """
        Instantiates the doppler shift model from redshift (z).

        Parameters
        -----------
        z : float
            the redshift
        """
        doppler_factor = z + 1
        return cls.from_doppler_factor(doppler_factor)

    @classmethod
    def from_doppler_factor(cls, doppler_factor):
        """
        Instantiates the doppler shift model from the doppler factor.

        Parameters
        -----------
        doppler_factor : float
            the doppler factor
        """
        velocity = constants.c*((doppler_factor**2 - 1)/(doppler_factor**2 + 1))
        return cls(velocity)

    def __init__(self, velocity):
        self.velocity = velocity

    def __call__(self, input):
        """
        Applies the doppler shift to the input, and returns the result

        Parameters
        -----------
        input : numpy array
            The input to be shifted
        """
        return self.__class__.evaluate(input, self.doppler_factor)

    @classmethod
    def evaluate(cls, input, doppler_factor):
        return input * doppler_factor

    def inverse(self):
        """
        Returns a new Doppler shift model, inverse of this doppler shift model.
        Basically, it instantiates the new doppler shift model with -velocity.
        """
        return DopplerShift(-self.velocity)

    @property
    def beta(self):
        return self.velocity / constants.c

    @property
    def doppler_factor(self):
        return math.sqrt((1 + self.beta)/(1 - self.beta))

    @property
    def redshift(self):
        return self.doppler_factor - 1


class MultispecIRAFCompositeWCS(BaseSpectrum1DWCS, CompositeWCS):
    """
    A specialized composite WCS model for IRAF FITS multispec specification.
    This model adds applies doppler adjustment and log (if applicable) to the
    dispersion WCS. The dispersion WCS may be a linear WCS or a weighted
    combination WCS. This model stores FITS metadata and also provides a method
    to generate the  multispec string for the FITS header.

    Parameters
    -----------
    dispersion_wcs : Model
        the wcs that returns the raw dispersion when passed the pixel
        coordinates. There is no doppler shift applied to this dispersion
    num_pixels : int
        the number of valid pixels
    z : float, optional
        the redshift, used to determine the doppler shift needed to correct the
        dispersion. Th default value indictes there is no shift
    log : bool, optional
        whether log correction needs to be applied to the dispersion
    aperture : int, optional
        the aperture number
    beam : int, optional
        the beam number
    aperture_low : float, optional
        the lower aperture limit
    aperture_high : float, optional
        the higher aperture limit
    unit: astropy.unit, optional
        the unit of the dispersion
    """

    def __init__(self, dispersion_wcs, num_pixels, z=1.0, log=False,
                 aperture=1, beam=88, aperture_low=0.0, aperture_high=0.0,
                 unit=None):
        doppler_wcs = DopplerShift.from_redshift(z).inverse
        super(MultispecIRAFCompositeWCS, self).__init__([dispersion_wcs,
                                                         doppler_wcs])
        self.aperture = aperture
        self.beam = beam
        self.num_pixels = num_pixels
        self.aperture_low = aperture_low
        self.aperture_high = aperture_high
        self.unit = unit
        if log:
            self.add_WCS(lambda x: 10 ** x)

    def __call__(self, pixel_indices):
        """
        Applies the model to the pixel coordinates, to produce the dispersion
        at those pixels

        Parameters
        -----------
        pixel_indices : numpy array
            the pixel coordinates on which dispersion needs to be computed
        """
        dispersion = super(MultispecIRAFCompositeWCS,
                           self).__call__(pixel_indices)
        return dispersion * self.unit

    def get_fits_spec(self):
        """
        Returns the list of spec parameters that represent this model in the
        order they are supposed to be stored in FITS files
        """
        if len(self.wcs_list) == 3:
            # we don't need to compute the log of dispersion while writing
            log_wcs = self.wcs_list.pop()
            dispersion_type = 1
        else:
            log_wcs = None
            dispersion_type = 0
        dispersion_wcs, doppler_wcs = self.wcs_list
        if isinstance(dispersion_wcs, WeightedCombinationWCS):
            dispersion_type = 2
        # else dispersion_wcs instance of Spectrum1DPolynomialWCS
        dispersion = self.__call__(np.arange(self.num_pixels)).value
        avg_disp_delta = (dispersion[1:] - dispersion[:-1]).mean()
        z = doppler_wcs.redshift
        spec = [self.aperture, self.beam, dispersion_type, dispersion[0],
                avg_disp_delta, self.num_pixels, z,
                self.aperture_low, self.aperture_high]
        if dispersion_type == 2:
            # add the parameters of the functions to the spec string
            for wcs, weight, zero_point_offset in dispersion_wcs.wcs_list:
                spec.extend([weight, zero_point_offset])
                spec.extend(wcs.get_fits_spec())
        if dispersion_type == 1:
            # add the log_wcs back
            self.add_WCS(log_wcs)

        return spec
