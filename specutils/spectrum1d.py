# Licensed under a 3-clause BSD style license - see LICENSE.rst
#This module implements the Spectrum1D class. It is eventually supposed to migrate to astropy core

from astropy.nddata import NDData

#!!!! checking scipy availability
scipy_available = True
try:
    import scipy
    from scipy import interpolate
except ImportError:
    scipy_available = False

import copy

import numpy as np

def spec_operation(func):
    def convert_operands(self, operand):
        if isinstance(operand, self.__class__):
            if not (all(self.dispersion == operand.dispersion) and\
                self.units == operand.units):
                raise ValueError('Dispersion and units need to match for both Spectrum1D objects')

        elif np.isscalar(operand):
            return func(self, operand)
        else:
            raise ValueError("unsupported operand type(s) for operation: %s and %s" %
                             (type(self), type(operand)))
    return convert_operands

class BoolMask(object):
    
    #Needs to live in NDData sooner or later
    #should also be a subclass of nddata, but am not sure how.
    
    def __init__(self, mask):
        self.mask = mask.astype(bool)
    
    def __getitem__(self, key):
        return BoolMask(self.mask[key])
    
    def interpolate():
        pass
    
    def bool_arithmetic(a, b):
        if isinstance(a, BoolMask) and isinstance(b, boolmask):
            return BoolMask(np.logical_and(a.mask, b.mask))
        else:
            raise ValueError("unsupported operand type(s): %s and %s" % \
                             (type(a), type(b)))
    
    def __add__(self, operand):
        return bool_arithmetic(self, operand)
    
    def __sub__(self, operand):
        return bool_arithmetic(self, operand)
        
    def __mul__(self, operand):
        return bool_arithmetic(self, operand)

    def __div__(self, operand):
        return bool_arithmetic(self, operand)
    
    
    #mirror functions
    def __radd__(self, operand):
        return self.__add__(operand, self)
        
    def __rsub__(self, operand):
        return self.__sub__(operand, self)
        
    def __rmul__(self, operand):
        return self.__mul__(operand, self)
            
    def __rdiv__(self, operand):
        return self.__div__(operand, self)


class Spectrum1D(NDData):
    """Class implementing a 1D spectrum"""
    
    def __init__(self, flux, dispersion=None, dispersion_unit=None,
                 error=None, mask=None, wcs=None, meta=None,
                 units=None, copy=True, validate=True):
        #needed to change order from (dispersion, flux) -> (flux, dispersion)
        #as dispersion=None for wcs.
        
        #added some WCS classes as I was not sure how to deal with both wcs and 
        
                
        
        NDData.__init__(self, data=flux, error=error, mask=mask,
                        wcs=wcs, meta=meta, units=units,
                        copy=copy, validate=validate)
        
        if wcs==None:
            self.dispersion = dispersion
            self.dispersion_unit = dispersion_unit
        else:
            self.wcs = wcs
            self.dispersion = wcs.get_lookup_table()
            self.dispersion_unit = wcs.units[0]
    
    @property
    def flux(self):
        #returning the flux
        return self.data
        
        
    #!!!! Not sure if we should have a setter for the flux
    #!!!! as we don't check if the new flux has the same shape as the error and mask
    
    #@flux.setter
    #def flux_setter(self, flux):
    #    self.data = flux
    
    
        
    def interpolate(self, dispersion, kind='linear', bounds_error=True, fill_value=np.nan, copy=True):
        """Interpolates onto a new wavelength grid and returns a `Spectrum1D`-obect
        Parameters:
        -----------
        
        dispersion: `numpy.ndarray`
            new dispersion array
        """
        
        if not scipy_available:
            if kind != 'linear':
                raise ValueError('Only \'linear\' interpolation is available if scipy is not installed')
            
            #### Erik & Thomas --- Can you tell me which logging stream to attach to and warn,
            #### about that bounds_error & fill_value is ignored as scipy not available
            interpolated_flux = np.interp(new_)
        else:
            spectrum_interp = interpolate.interp1d(self.disp, self.flux,
                                        kind=kind, bounds_error=bounds_error,
                                        fill_value=fill_value)
            new_flux = spectrum_interp(dispersion)
            
            if error!=None:
                new_error = self.error.interpolate(dispersion,
                                                   kind=kind,
                                                   bounds_error=bounds_error,
                                                   fill_value=fill_value)
            else:
                new_error = None
            
            if mask!=None:
                new_mask = self.mask.interpolate(dispersion,
                                                   kind=kind,
                                                   bounds_error=bounds_error,
                                                   fill_value=fill_value)
            else:
                new_mask = None
            
            
            
        
        if copy:    
            return self.__class__(new_flux, dispersion, error=new_error,
                                mask=new_mask, meta=copy.deepcopy(meta),
                                copy=False)
        else:
            raise NotImplementedError('Inplace will be implemented soon')
        
        
        
    def slice(self, start=None, stop=None, units='dispersion'):
        """Slicing the spectrum
        
        Paramaters:
        -----------
        
        start: numpy.float or int
            start of slice
        stop:  numpy.float or int
            stop of slice
        units: str
            allowed values are 'disp', 'pixel'
        """
        
        if units == 'dispersion':
            if start == None:
                start_idx = None
            else:
                start_idx = self.dispersion.searchsorted(start)
            
            if stop == None:
                stop_idx = None
            else:
                stop_idx = self.disp.searchsorted(stop)
            
            spectrum_slice = slice(start_idx, stop_idx)
        elif units == 'index':
            spectrum_slice = slice(start, stop)
        else:
            raise ValueError("units keyword can only have the values 'disp', 'index'")
        
        
        if copy:
            return self.__class__(self.flux[spectrum_slice],
                                  self.dispersion[spectrum_slice],
                                  error=self.error[spectrum_slice],
                                  mask=self.mask[spectrum_slice],
                                  meta=copy.deepcopy(meta))
        else:
            raise NotImplementedError('Inplace will be implemented soon')


