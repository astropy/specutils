# Licensed under a 3-clause BSD style license - see LICENSE.rst
#This module implements the Spectrum1D class. It is eventually supposed to migrate to astropy core

from astropy.nddata import NDData
from scipy import interpolate
import numpy as np

def spec_operation(func):
    def convert_operands(self, operand):
        if isinstance(operand, self.__class__):
            if all(self.disp == operand.disp):
                return func(self, operand.flux)
            else:
                new_disp = np.union1d(self.disp, operand.disp)
                return func(self.interpolate(new_disp), operand.interpolate(new_disp).flux)

        elif np.isscalar(operand):
            return func(self, operand)
        else:
            raise ValueError("unsupported operand type(s) for operation: %s and %s" %
                             (type(self), type(operand)))
    return convert_operands


class Spectrum1D(NDData):
    """Class implementing a 1D spectrum"""
    
    @classmethod
    def from_dispflux(cls, disp, flux, error=None, mask=None):
        """Initializing `Spectrum1D`-object from two `numpy.ndarray` object
        
        Paramateres:
        ------------
            disp: `numpy.ndarray`
                dispersion solution (e.g. wavelength array)
            flux: `numpy.ndarray`
                flux array"""
        if disp.ndim != 1 or disp.shape != flux.shape:
            raise ValueError('disp and flux need to be one-dimensional arrays with the same shape')
            
        return cls(data=flux, wcs=disp, error=error, mask=mask)
    
    @classmethod
    def from_fits(cls, filename, error=None):
        """This is an example function to demonstrate how
        classmethods are a clean way to instantiate Spectrum1D objects"""
        raise NotImplementedError('This function is not implemented yet')
    
    
    @property
    def flux(self):
        #returning the flux
        return self.data
        
    @flux.setter
    def flux_setter(self, flux):
        self.data = flux
    
    @property
    def disp(self):
        #returning the disp
        return self.wcs
    
        
    def interpolate(self, new_disp, kind='linear', bounds_error=True, fill_value=np.nan):
        """Interpolates onto a new wavelength grid and returns a `Spectrum1D`-obect
        Parameters:
        -----------
        
        new_disp: `numpy.ndarray`
            new dispersion array
        """
        spectrum_interp = interpolate.interp1d(self.disp, self.flux,
                                        kind=kind, bounds_error=bounds_error,
                                        fill_value=fill_value)
        new_flux = spectrum_interp(new_disp)
        
        return self.__class__.from_dispflux(new_disp, new_flux)
        
        
    def slice(self, start=None, stop=None, units='disp'):
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
        
        if units == 'disp':
            if start == None:
                start_idx = None
            else:
                start_idx = self.disp.searchsorted(start)
            
            if stop == None:
                stop_idx = None
            else:
                stop_idx = self.disp.searchsorted(stop)
            
            spectrum_slice = slice(start_idx, stop_idx)
        elif units == 'pixel':
            spectrum_slice = slice(start, stop)
        else:
            raise ValueError("units keyword can only have the values 'disp', 'pixel'")
        
        return self.__class__.from_dispflux(self.disp[spectrum_slice], self.flux[spectrum_slice])
    
    @spec_operation
    def __add__(self, operand):
        
        """Adds two spectra together, or adds finite real numbers across an entire spectrum."""
        
        return self.__class__.from_dispflux(self.disp, self.flux + operand)
        

    @spec_operation
    def __sub__(self, operand):
        
        """Subtracts two spectra, or subtracts a finite real numbers from an entire spectrum."""
        
        return self.__class__.from_dispflux(self.disp, self.flux - operand)
        

    @spec_operation
    def __mul__(self, operand):
        
        """Multiplies two spectra, or multiplies a finite real numbers across an entire spectrum."""
        
        return self.__class__.from_dispflux(self.disp, self.flux * operand)
        

    @spec_operation
    def __div__(self, operand):
        
        """Divides two spectra, or divides a finite real numbers across an entire spectrum."""
        
        return self.__class__.from_dispflux(self.disp, self.flux / operand)
        
    @spec_operation
    def __pow__(self, operand):
        
        """Performs power operations on spectra."""
        
        return self.__class__.from_dispflux(self.disp, self.flux ** operand)
        

    def __len__(self):
        return len(self.disp)


    # Mirror functions
    
    def __radd__(self, spectrum, **kwargs):
        return self.__add__(spectrum, **kwargs)
        
    def __rsub__(self, spectrum, **kwargs):
        return self.__sub__(spectrum, **kwargs)
        
    def __rmul__(self, spectrum, **kwargs):
        return self.__mul__(spectrum, **kwargs)
            
    def __rdiv__(self, spectrum, **kwargs):
        return self.__div__(spectrum, **kwargs)
    
    def __rpow__(self, spectrum, **kwargs):
        return self.__pow__(spectrum, **kwargs)
        