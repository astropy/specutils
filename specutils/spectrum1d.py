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
        elif np.isscalar(operand):
            flux = operand
            
        else:
            raise ValueError("unsupported operand type(s) for operation: %s and %s" %
                             (type(self), type(operand)))
        
        return func(self, flux)
        
    return convert_operands



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
    def __add__(self, operand_flux):

        new_flux = self.flux + operand_flux
        new_meta = merge_meta(self.meta, operand_meta)
        return self.__class__(new_flux,
                              #!!! What if it's a WCS
                              self.dispersion.copy(),
                              meta=new_meta)
