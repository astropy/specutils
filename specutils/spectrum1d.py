# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module implements the Spectrum1D class.

from __future__ import print_function, division

__all__ = ['Spectrum1D']

from astropy import log
from astropy.nddata import NDData, FlagCollection

from astropy.utils import misc

from specutils.wcs import BaseSpectrum1DWCS, Spectrum1DLookupWCS, Spectrum1DLinearWCS, Spectrum1DWCSFITSError


from astropy.io import fits
from astropy import units as u

import numpy as np


class Spectrum1D(NDData):
    """A subclass of `NDData` for a one dimensional spectrum in Astropy.
    
    This class inherits all the base class functionality from the NDData class
    and is communicative with other Spectrum1D objects in ways which make sense.

    Parameters
    ----------
    data: array
    wcs: wcs
    unit: unit
    """

    _wcs_attributes = {'wavelength': {'unit': u.m},
                      'frequency': {'unit': u.Hz},
                      'energy': {'unit': u.J}}
    
    @classmethod
    def from_array(cls, dispersion, flux, dispersion_unit=None, uncertainty=None, mask=None,
                   flags=None, meta=None, copy=True,
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

        flags : `~numpy.ndarray` or `~astropy.nddata.FlagCollection`, optional
            Flags giving information about each pixel. These can be specified
            either as a Numpy array of any type with a shape matching that of the
            data, or as a `~astropy.nddata.FlagCollection` instance which has a
            shape matching that of the data.

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
            raise ValueError("dispersion and flux need to be one-dimensional Numpy arrays with the same shape")
        spec_wcs = Spectrum1DLookupWCS(dispersion, unit=dispersion_unit)

        if copy:
            flux = flux.copy()

        return cls(data=flux, wcs=spec_wcs, unit=unit, uncertainty=uncertainty,
                   mask=mask, flags=flags, meta=meta)
    
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

        if isinstance(flag_columns, basestring):
            flags = table[flag_columns]
        elif misc.isiterable(flag_columns):
            flags = FlagCollection(shape=flux.shape)
            for flag_column in flag_columns:
                flags[flag_column] = table[flag_column]
        else:
            raise ValueError('flag_columns should either be a string or a list (or iterable) of strings')

        return cls.from_array(flux=flux.data, dispersion=dispersion.data,
                              uncertainty=uncertainty, dispersion_unit=dispersion.units,
                              unit=flux.units, mask=table.mask, flags=flags,
                              meta=table.meta)
        
    
    
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
    def from_fits(cls, filename, uncertainty=None):
        """This is an example function to demonstrate how
        classmethods are a clean way to instantiate Spectrum1D objects"""
        header = fits.getheader(filename)
        try:
            cls.dispersion = Spectrum1DLinearWCS.from_header(header)
        except:
            pass
        raise NotImplementedError('This function is not implemented yet')

    def __init__(self, *args, **kwargs):
        super(Spectrum1D, self).__init__(*args, **kwargs)
        for key in self._wcs_attributes:

            wcs_attribute_unit = self._wcs_attributes[key]['unit']
            if wcs_attribute_unit.physical_type == self.wcs.unit.physical_type:
                wcs_attribute_unit = self.wcs.unit

            if not wcs_attribute_unit.is_equivalent_to(self.wcs.unit):
                del self._wcs_attributes[key]

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
        return self.__dict__.keys() + self._wcs_attributes.keys() + \
               [item + '_unit' for item in self._wcs_attributes.keys()]

    
    @property
    def flux(self):
        #returning the flux
        return self.data
        
    @flux.setter
    def flux_setter(self, flux):
        self.data = flux
    
    @property
    def dispersion(self):
        #returning the disp
        if not hasattr(self.wcs, 'lookup_table'):
            self.wcs.lookup_table = self.wcs(np.arange(len(self.flux)))

        return self.wcs.lookup_table

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
        
        # Transform the dispersion end points to index space
        start_index, stop_index = self.wcs([start, stop])
        
        #return self.slice_index(start_index, stop_index)
    
    
    def slice_index(self, start=None, stop=None):
        """Slice the spectrum within a given start and end index.
        
        Parameters
        ----------
        start : `float`
            Starting slice point.
        stop : `float`
            Stopping slice point.
        
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
        # >> disp, flux, error, mask, and flags
        # Which are all common NDData objects, therefore I am (perhaps
        # reasonably) assuming that __slice__ will be a NDData base function
        # which we will inherit.
        raise NotImplementedError('Will presumeably implemented in core NDDATA,'
                                  'though this is just trivial indexing.')
        return self[start:stop]

    def to_fits(self, filename):
        """
        Write to fits file

        Parameters
        ----------

        filename : `str`
            file name to write the current spectrum to in FITS format

        Notes
        -----


        """
        fits_file = fits.PrimaryHDU(self.data)
        try:
            self.wcs.to_fits_header(fits_file.header)
        except AttributeError:
            raise Spectrum1DWCSFITSError('Current WCS does not support writing to FITS files - interpolate to one that '
                                         'does')

        fits_file.writeto(filename)
