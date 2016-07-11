# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module implements the Spectrum1D class.

from __future__ import print_function, division
from specutils.models.Indexer import Indexer

__all__ = ['Spectrum1D']

import copy
from astropy.extern import six
from astropy import log
from astropy.nddata import NDData, FlagCollection

from astropy.utils import misc

from specutils.wcs import BaseSpectrum1DWCS, Spectrum1DLookupWCS


from astropy import units as u

import numpy as np


class Spectrum1D(NDData):
    """A subclass of `NDData` for a one dimensional spectrum in Astropy.

    This class inherits all the base class functionality from the NDData class
    and is communicative with other Spectrum1D objects in ways which make sense.

    Parameters
    ----------
    data : `~numpy.ndarray`
        flux of the spectrum

    wcs : `spectrum1d.wcs.specwcs.BaseSpectrum1DWCS`-subclass
        transformation between pixel coordinates and "dispersion" coordinates
        this carries the unit of the dispersion

    unit : `~astropy.unit.Unit` or None, optional
        unit of the flux, default=None

    mask : `~numpy.ndarray`, optional
        Mask for the data, given as a boolean Numpy array with a shape
        matching that of the data. The values must be ``False`` where
        the data is *valid* and ``True`` when it is not (like Numpy
        masked arrays). If `data` is a numpy masked array, providing
        `mask` here will causes the mask from the masked array to be
        ignored.

    meta : `dict`-like object, optional
        Metadata for this object.  "Metadata" here means all information that
        is included with this object but not part of any other attribute
        of this particular object.  e.g., creation date, unique identifier,
        simulation parameters, exposure time, telescope name, etc.

    """

    _wcs_attributes = {'wavelength': {'unit': u.m},
                      'frequency': {'unit': u.Hz},
                      'energy': {'unit': u.J},
                      'velocity': {'unit': u.m/u.s}}

    @classmethod
    def from_array(cls, dispersion, flux, dispersion_unit=None,
                   uncertainty=None, mask=None, meta=None, copy=True,
                   unit=None):
        """Initialize `Spectrum1D`-object from two `numpy.ndarray` objects

        Parameters:
        -----------
        dispersion : `~astropy.units.quantity.Quantity` or `~np.array`
            The dispersion for the Spectrum (e.g. an array of wavelength
            points). If an array is specified `dispersion_unit` needs to be a spectral unit

        flux : `~astropy.units.quantity.Quantity` or `~np.array`
            The flux level for each wavelength point. Should have the same length
            as `dispersion`.

        dispersion_unit :
        error : `~astropy.nddata.NDError`, optional
            Errors on the data.

        mask : `~numpy.ndarray`, optional
            Mask for the data, given as a boolean Numpy array with a shape
            matching that of the data. The values should be ``False`` where the
            data is *valid* and ``True`` when it is not (as for Numpy masked
            arrays).

        meta : `dict`-like object, optional
            Metadata for this object. "Metadata here means all information that
            is included with this object but not part of any other attribute
            of this particular object. e.g., creation date, unique identifier,
            simulation parameters, exposure time, telescope name, etc.

        copy : bool, optional
            If True, the array will be *copied* from the provided `data`,
            otherwise it will be referenced if possible (see `numpy.array` :attr:`copy`
            argument for details).

        Raises
        ------
        ValueError
            If the `dispersion` and `flux` arrays cannot be broadcast (e.g. their shapes
            do not match), or the input arrays are not one dimensional.

        """

        if dispersion.ndim != 1 or dispersion.shape != flux.shape:
            raise ValueError("dispersion and flux need to be one-dimensional "
                             "Numpy arrays with the same shape")

        if hasattr(dispersion, 'unit'):
            if dispersion_unit is not None:
                dispersion = dispersion.to(dispersion_unit).value
            else:
                dispersion_unit = dispersion.unit
                dispersion = dispersion.value


        spec_wcs = Spectrum1DLookupWCS(dispersion, unit=dispersion_unit)

        if copy:
            flux = flux.copy()

        return cls(flux=flux, wcs=spec_wcs, unit=unit, uncertainty=uncertainty,
                   mask=mask, meta=meta)

    @classmethod
    def from_table(cls, table, dispersion_column='dispersion',
                   flux_column='flux', uncertainty_column=None,
                   flag_columns=None):
        """
        Initializes a `Spectrum1D`-object from an `~astropy.table.Table` object

        Parameters
        ----------

        table : ~astropy.table.Table object

        dispersion_column : str, optional
            name of the dispersion column. default is 'dispersion'

        flux_column : str, optional
            name of the flux column. default is 'flux'

        uncertainty_column : str, optional
            name of the uncertainty column. If set to None uncertainty is set to None. default is None

        flag_columns : str or list, optional
            name or names of flag columns. If multiple names are supplied a ~astropy.nddata.FlagCollection will be built.
            default is None
        """

        flux = table[flux_column]
        dispersion = table[dispersion_column]

        if uncertainty_column is not None:
            uncertainty = table[uncertainty_column]
            if uncertainty.unit != flux.unit:
                log.warning('"uncertainty"-column and "flux"-column do not share the units (%s vs %s) ',
                            uncertainty.unit, flux.unit)
        else:
            uncertainty = None

        return cls.from_array(flux=flux.data, dispersion=dispersion.data,
                              uncertainty=uncertainty, dispersion_unit=dispersion.units,
                              unit=flux.units, mask=table.mask, meta=table.meta)



    @classmethod
    def from_ascii(cls, filename, uncertainty=None, mask=None, dtype=np.float, comments='#',
                   delimiter=None, converters=None, skiprows=0,
                   usecols=None):
        raw_data = np.loadtxt(filename, dtype=dtype, comments=comments,
                              delimiter=delimiter, converters=converters,
                              skiprows=skiprows, usecols=usecols, ndmin=2)

        if raw_data.shape[1] != 2:
            raise ValueError('data contained in filename must have exactly two columns')

        return cls.from_array(dispersion=raw_data[:,0], flux=raw_data[:,1], uncertainty=uncertainty, mask=mask)

    @classmethod
    def from_fits(cls, filename):
        """
        This function is a dummy function and will fail for now. Please use the functions provided in
        `~specutils.io.read_fits` for this task.
        """

        raise NotImplementedError('This function is not implemented. To read FITS files please refer to the'
                                  ' documentation')

    def __init__(self, flux, wcs, unit=None, uncertainty=None, mask=None,
                 meta=None, indexer=None, *args, **kwargs):

        super(Spectrum1D, self).__init__(data=flux, unit=unit, wcs=wcs, uncertainty=uncertainty,
                   mask=mask, meta=meta, *args, **kwargs)

        self._wcs_attributes = copy.deepcopy(self.__class__._wcs_attributes)
        if indexer is None:
            self.indexer = Indexer(0, len(flux))
        else:
            self.indexer = indexer
        for key in list(self._wcs_attributes):

            wcs_attribute_unit = self._wcs_attributes[key]['unit']

            try:
                unit_equivalent = wcs_attribute_unit.is_equivalent(self.wcs.unit, equivalencies=self.wcs.equivalencies)
            except TypeError:
                unit_equivalent = False


            if not unit_equivalent:
                #if unit is not convertible to wcs attribute - delete that wcs attribute
                del self._wcs_attributes[key]
                continue

            if wcs_attribute_unit.physical_type == self.wcs.unit.physical_type:
                self._wcs_attributes[key]['unit'] = self.wcs.unit


    def flux_getter(self):
        #returning the flux
        return u.Quantity(self.data, self.unit, copy=False)

    def flux_setter(self, flux):
        if hasattr(flux, 'unit'):
            if self.unit is not None:
                flux = flux.to(self.unit).value
            else:
                raise ValueError('Attempting to set a new unit for this object'
                                 'this is not allowed by Spectrum1D')

        self._data = flux


    flux = property(flux_getter, flux_setter)

    def __getattr__(self, name):
        if name in self._wcs_attributes:
            return self.dispersion.to(self._wcs_attributes[name]['unit'], equivalencies=self.wcs.equivalencies)
        elif name[:-5] in self._wcs_attributes and name[-5:] == '_unit':
            return self._wcs_attributes[name[:-5]]['unit']
        else:
            super(Spectrum1D, self).__getattribute__(name)


    def __setattr__(self, name, value):
        if name[:-5] in self._wcs_attributes and name[-5:] == '_unit':
            self._wcs_attributes[name[:-5]]['unit'] = u.Unit(value)
        else:
            super(Spectrum1D, self).__setattr__(name, value)

    def __dir__(self):
        return list(self.__dict__.keys()) + list(self._wcs_attributes.keys()) + \
               [item + '_unit' for item in self._wcs_attributes.keys()]




    #TODO: let the WCS handle what to do with len(flux)
    @property
    def dispersion(self):
        #returning the disp
        pixel_indices = np.arange(len(self.flux))
        return self.wcs(self.indexer(pixel_indices))

    @property
    def dispersion_unit(self):
        return self.wcs.unit

        
    def interpolate(self, new_dispersion, kind='linear', bounds_error=True, fill_value=np.nan):
        """Interpolates onto a new wavelength grid and returns a new `Spectrum1D`-object.
        
        Parameters
        ----------
        new_dispersion : `~numpy.ndarray`
            The dispersion array to interpolate the flux on to.
        
        kind : `str` or `int`, optional
            Specifies the kind of interpolation as a string
            ('linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic')
            or as an integer specifying the order of the spline interpolator
            to use. Default is 'linear'.
        
        bounds_error : `bool`, optional
            If True, an error is thrown any time interpolation is attempted on a
            dispersion point outside of the range of the original dispersion map
            (where extrapolation is necessary). If False, out of bounds values
            are assigned `fill_value`. By default, an error is raised.
            
        fill_value : `float`, optional
            If provided, then this value will be used to fill in for requested
            dispersion points outside of the original dispersion map. If not
            provided, then the default is NaN.
        
        Raises
        ------
        ImportError
            If the `SciPy interpolate interp1d <http://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html>`_
            function cannot be imported.
            
        Notes
        -----
        When the `Spectrum1D` class has an associated error array, the nearest
        uncertainty is taken for each new dispersion point.
        
        """
        
        # Check for SciPy availability


        if kind != 'linear':
            raise ValueError('No other kind but linear supported')

        if not isinstance(new_dispersion, BaseSpectrum1DWCS):
            new_dispersion = Spectrum1DLookupWCS(np.array(new_dispersion))


        new_pixel = self.wcs.invert(new_dispersion.lookup_table)

        new_flux = np.interp(new_pixel, self.wcs.pixel_index, self.flux, left=np.nan, right=np.nan)
        

        return self.__class__(new_flux, wcs=new_dispersion, meta=self.meta)

        
    def slice_dispersion(self, start=None, stop=None):
        """Slice the spectrum within a given start and end dispersion value.
        
        Parameters
        ----------
        start : `float`
            Starting slice point.
        stop : `float`
            Stopping slice point.
        
        Notes
        -----
        Often it is useful to slice out a portion of a `Spectrum1D` objects
        either by two dispersion points (e.g. two wavelengths) or by the indices
        of the dispersion/flux arrays (see :meth:`~Spectrum1D.slice_index` for this
        functionality).
        
        Examples
        --------
        
        >>> from specutils import Spectrum1D
        >>> from astropy import units
        >>> import numpy as np
        >>> dispersion = np.arange(4000, 5000, 0.12)
        >>> flux = np.random.randn(len(dispersion))
        >>> mySpectrum = Spectrum1D.from_array(dispersion,
                                               flux,
                                               dispersion_unit=units.m)
        
        >>> # Now say we wanted a slice near H-beta at 4861 Angstroms
        >>> hBeta = mySpectrum.slice_dispersion(4851.0, 4871.0)
        >>> hBeta
        <hBeta __repr__ #TODO>
        
        See Also
        --------
        See `~Spectrum1D.slice_index`
        """
        raise NotImplementedError('Waiting for slicing implementation in WCS and NDData')
        # Transform the dispersion end points to index space
        start_index, stop_index = self.wcs([start, stop])
        
        #return self.slice_index(start_index, stop_index)
    
    
    def slice_index(self, start=None, stop=None, step=None):
        """Slice the spectrum within a given start and end index.
        
        Parameters
        ----------
        start : int
            Starting slice point.
        stop : int
            Stopping slice point.
        step : int
            Slice step
        Notes
        -----
        Often it is useful to slice out a portion of a `Spectrum1D` objects
        either by two index points (see :meth:`~Spectrum1D.slice_dispersion`) or by
        the indices of the dispersion/flux array.
        
        See Also
        --------
        See `~Spectrum1D.slice_dispersion`
        """
        
        # We need to slice the following items:
        # >> disp, flux, error, and mask
        # Which are all common NDData objects, therefore I am (perhaps
        # reasonably) assuming that __slice__ will be a NDData base function
        # which we will inherit.
        # At this time, that function raises an error if WCS is not None, so it
        # cannot be used
        item = slice(start, stop, step)
        new_data = self.data[item]

        if self.uncertainty is not None:
            new_uncertainty = self.uncertainty[item]
        else:
            new_uncertainty = None

        if self.mask is not None:
            new_mask = self.mask[item]
            # mask setter expects an array, always
            if new_mask.shape == ():
                new_mask = np.array(new_mask)
        else:
            new_mask = None

        new_indexer = self.indexer.__getitem__(item)
        new_wcs = self.wcs

        return self.__class__(new_data, new_wcs, meta=self.meta, unit=self.unit
                              , uncertainty=new_uncertainty, mask=new_mask,
                              indexer=new_indexer)


