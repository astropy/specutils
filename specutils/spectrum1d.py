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


def warn2logging(message, logstream=None):
    #!!!!!! Thomas you are the logging expert, I will call this function whenever I want to log to warn.
    #!!!!! Can you tell me what to write in here.
    print 'WARNING: %s' % message
    
    
    
def merge_meta(meta1, meta2, logstream=None):
    #Merging meta information and removing all duplicate keys -- warning what keys were removed
    #should be in NDData somewhere
    meta1_keys = meta1.viewkeys()
    meta2_keys = meta2.viewkeys()
    
    
    duplicates = meta1_keys & meta2_keys
    
    if len(duplicates) > 0:
        warn2logging('Removing duplicate keys found in meta data: ' + ','.join(duplicates), logstream=logstream)
    
    new_meta = copy.deepcopy(meta1)
    new_meta.update(copy.deepcopy(meta2))
    for key in duplicates:
        del new_meta[key]
    
    return new_meta

def spec_operation(func):
    #used as a decorator for the arithmetic of spectra
    def convert_operands(self, operand):
        
        #checking if they have the same wcs and units
        if isinstance(operand, self.__class__):
            if not (all(self.dispersion == operand.dispersion) and\
                self.units == operand.units):
                raise ValueError('Dispersion and units need to match for both Spectrum1D objects')
            
            flux = operand.flux
            mask = operand.mask
            error = operand.error
            meta = operand.meta
        #for a scalar the flux is the scalar and the error and mask and meta are None/ {}
        elif np.isscalar(operand):
            flux = operand
            mask = None
            error = None
            meta = {}
            
        else:
            raise ValueError("unsupported operand type(s) for operation: %s and %s" %
                             (type(self), type(operand)))
        
        return func(self, flux, error, mask, meta)
        
    return convert_operands


class BoolMask(np.ndarray):
    # !!! To be discussed; should probably live somewhere near NDData
    def __new__(cls, input_array):
        obj = np.asarray(input_array, dtype='bool').view(cls)
        return obj
    
    def interpolate(self, old_lookup_table, new_lookup_table):
        new_mask_raw = np.interp(new_lookup_table, old_lookup_table, self.astype(float64), left=1, right=1)
        return np.ceil(new_mask_raw).astype(bool)
    
    def bool_arithmetic(self, operand):
        # !!!! What happens if None?
        if operand == None:
            return self
        if not isinstance(operand, BoolMask):
            raise ValueError('unsupported operand type(s) for +: %s and %s' %\
                             (type(self), type(operand)))
        return np.logical_or(self, operand)

    def mask_add(self, operand):
        return self.bool_arithmetic(operand)
        
    def mask_sub(self, operand):
        return self.bool_arithmetic(operand)
    
    def mask_mul(self, operand):
        return self.bool_arithmetic(operand)
        
    def mask_div(self, operand):
        return self.bool_arithmetic(operand)


class SDError(np.ndarray):
    
    def __new__(cls, input_array):
        obj = np.asarray(input_array).view(cls)
        return obj
    
    def interpolate(self, old_lookup_table, new_lookup_table):
        return np.interp(new_lookup_table, old_lookup_table)
    
    def check_operand(self, operand):
    #checking if both of the operands are SDerrors or raising an exception
    
    # !!!! What happens if None?
        if not isinstance(operand, SDError):
            raise ValueError('unsupported operand type(s) for +: %s and %s' %\
                             (type(self), type(operand)))
    
    def error_add(self, self_data, operand, operand_data, result_data):
        self.check_operand(operand)
        return np.sqrt(self**2 + operand**2)
        
    def error_sub(self, self_data, operand, operand_data, result_data):
        self.check_operand(operand)
        return np.sqrt(self**2 + operand**2)

    def error_mul(self, self_data, operand, operand_data, result_data):
        self.check_operand(operand)
        return np.sqrt((self / self_data)**2 + (operand / operand_data)**2)
        
    def error_div(self, self_data, operand, operand_data, result_data):
        self.check_operand(operand)
        return np.sqrt((self / self_data)**2 + (operand / operand_data)**2)



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
            spectrum_interp = interpolate.interp1d(self.dispersion, self.flux,
                                        kind=kind, bounds_error=bounds_error,
                                        fill_value=fill_value)
            new_flux = spectrum_interp(dispersion)
            
            if error!=None:
                new_error = self.error.interpolate(self.dispersion, dispersion,
                                                   kind=kind,
                                                   bounds_error=bounds_error,
                                                   fill_value=fill_value)
            else:
                new_error = None
            
            if mask!=None:
                new_mask = self.mask.interpolate(self.dispersion, dispersion,
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
        
        
    #naming convention start and stop taken from python slices. these are nothing else but slices
    #so i think we should keep start stop step
    def slice(self, start=None, stop=None, step=None, units='dispersion'):
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
            if step != None:
                raise ValueError('step can only be specified for units=pixel')
            if start == None:
                start_idx = None
            else:
                start_idx = self.dispersion.searchsorted(start)
            
            if stop == None:
                stop_idx = None
            else:
                stop_idx = self.dispersion.searchsorted(stop)
            
            spectrum_slice = slice(start_idx, stop_idx)
        elif units == 'index':
            spectrum_slice = slice(start, stop, step)
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

    @spec_operation
    def __add__(self, operand_flux, operand_error, operand_mask, operand_meta):
        new_flux = self.flux + operand_flux
        
        if self.error != None and operand_error!=None:
            new_error = self.error.error_add(self.flux, operand_error, operand_flux, new_flux)
        else:
            new_error = None
        
        if self.mask != None:
            new_mask = self.mask.mask_add(operand_mask)
        elif operand_mask != None:
            new_mask = operand_mask.copy()
        else:
            new_mask = None
        
        new_meta = merge_meta(self.meta, operand_meta)
        
        return self.__class__(new_flux,
                              #!!! What if it's a WCS
                              self.dispersion.copy(),
                              error=new_error,
                              mask=new_mask,
                              meta=new_meta)
    
    @spec_operation
    def __sub__(self, operand_flux, operand_error, operand_mask, operand_meta):
        new_flux = self.flux + operand_flux
        
        if self.error != None and operand_error!=None:
            new_error = self.error.error_sub(self.flux, operand_error, operand_flux, new_flux)
        else:
            new_error = None
        
        if self.mask != None:
            new_mask = self.mask.mask_sub(operand_mask)
        elif operand_mask != None:
            new_mask = operand_mask.copy()
        else:
            new_mask = None
        
        new_meta = merge_meta(self.meta, operand_meta)
        
        return self.__class__(new_flux,
                              #!!! What if it's a WCS
                              self.dispersion.copy(),
                              error=new_error,
                              mask=new_mask,
                              meta=new_meta)
    
    @spec_operation
    def __mul__(self, operand_flux, operand_error, operand_mask, operand_meta):
        new_flux = self.flux + operand_flux
        
        if self.error != None and operand_error!=None:
            new_error = self.error.error_mul(self.flux, operand_error, operand_flux, new_flux)
        else:
            new_error = None
        
        if self.mask != None:
            new_mask = self.mask.mask_mul(operand_mask)
            
        elif operand_mask != None:
            new_mask = operand_mask.copy()
        else:
            new_mask = None
        
        new_meta = merge_meta(self.meta, operand_meta)
        
        return self.__class__(new_flux,
                              #!!! What if it's a WCS
                              self.dispersion.copy(),
                              error=new_error,
                              mask=new_mask,
                              meta=new_meta)
    
    @spec_operation
    def __div__(self, operand_flux, operand_error, operand_mask, operand_meta):
        new_flux = self.flux + operand_flux
        
        if self.error != None and operand_error!=None:
            new_error = self.error.error_div(self.flux, operand_error, operand_flux, new_flux)
        else:
            new_error = None
        
        if self.mask != None:
            new_mask = self.mask.mask_div(operand_mask)
        elif operand_mask != None:
            new_mask = operand_mask.copy()
        else:
            new_mask = None
        
        new_meta = merge_meta(self.meta, operand_meta)
        
        return self.__class__(new_flux,
                              #!!! What if it's a WCS
                              self.dispersion.copy(),
                              error=new_error,
                              mask=new_mask,
                              meta=new_meta)

if __name__ == '__main__':
    my_first_spec = Spectrum1D(np.random.rand(1000),
                               linspace(4000, 7000, 1000),
                               error=SDError(np.random.rand(1000)),
                               mask = BoolMask(np.random.rand(1000) > .5),
                               meta = dict(a=5, b=7, c=9))
    
    my_second_spec = Spectrum1D(np.random.rand(1000),
                               linspace(4000, 7000, 1000),
                               error=SDError(np.random.rand(1000)),
                               mask = BoolMask(np.random.rand(1000) > .5),
                               meta = dict(c=9, d=11, e=13))

