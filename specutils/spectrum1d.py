# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module implements the Spectrum1D class. It is eventually supposed to migrate to astropy core

from astropy.nddata import NDData
from scipy import interpolate
import numpy as np  

class SpectrumErrorArray(object):
    """ Represents the error axis of a Spectrum object
    
    .. warning::
        This is skeleton code!
    
    Spectral data errors can be represented in a variety of ways, e.g. 1-sigma
    errors, inverse variance, variance, etc. In order to provide generic functionality,
    the error data must carry with it some metadata and operations that help
    handle transformations of the error array.

    Parameters
    ----------
    array : `~numpy.ndarray`
        The error values
    
    Notes
    -----
    This class should be 'private' in that it should really only be subclassed
    and used by specific classes that implement

    """
    def __init__(self, values):
        self.values = values
    
    @abstractmethod
    def interpolate(self, new_dispersion):
        pass

class InverseVarianceErrorArray(SpectrumErrorArray):
    # ...
    
    def interpolate(self, new_dispersion):
        pass

class VarianceErrorArray(SpectrumErrorArray):
    # ...
    
    def interpolate(self, new_dispersion):
        pass

class Spectrum1DBase(NDData):
    """ Base class for Spectrum1D objects.
    
    .. warning::
        This is skeleton code!
    
    `Spectrum1DBase` is a stripped-down superclass of the Spectrum1D object
    that allows for creation of a Spectrum-like object **without** an internal
    dispersion array. This allows for the possibility of creating a 
    SpectrumCollection type class that could contain many spectra that share
    the same dispersion axis (e.g. SDSS spectra from the same plate).

    Parameters
    ----------
    flux : `~numpy.ndarray`
        The flux data as an array
        
    Notes
    -----
    This class should be 'private' in that it should really only be subclassed
    and used by Spectrum1D and SpectrumCollection

    """
    
    pass
    
class Spectrum1D(Spectrum1DBase):
    """ Class for 1-dimensional Spectrum objects.
    
    .. warning::
        This is skeleton code!
    
    `Spectrum1D` provides a container for 1-dimensional spectral data as well
    as generic operations for this data.

    Parameters
    ----------
    flux : `~numpy.ndarray`
        The flux data as an array
    dispersion : `~numpy.ndarray`
        An array of data representing the dispersion axis of the spectral data,
        e.g. wavelength, frequency, energy, velocity, etc.
    <NDData parameters...>
        
    Notes
    -----
    
    """
    
    
    def slice_dispersion(self, start=None, end=None):
        """ Slices the data arrays based on values in the dispersion array
        
            Parameters
            ----------
            start : any numeric type, optional
                The minimum value in the dispersion axis to slice on
            end : any numeric type, optional
                The maximum value in the dispersion axis to slice on
        """
        pass
    
    def slice_pixel(self, start=None, end=None):
        """ Slices the data arrays on pixel index (e.g. array index)
        
            Parameters
            ----------
            start : int, optional
                The starting index to slice the arrays on
            end : int, optional
                The ending index to slice the arrays on
            
            Notes
            -----
            This is equivalent to slicing each array like array[start:end]
            
        """
        pass
    
    def interpolate_linear(self, new_dispersion):
        """ Uses linear interpolation to resample the internal arrays onto
            the specified dispersion axis.
        
            Parameters
            ----------
            new_dispersion : `~numpy.ndarray`
                The new dispersion array to interpolate on to
            
        """
        pass
        
    def interpolate_spline(self, new_dispersion):
        """ Uses spline interpolation to resample the internal arrays onto
            the specified dispersion axis.
        
            Parameters
            ----------
            new_dispersion : `~numpy.ndarray`
                The new dispersion array to interpolate on to
            
        """
        pass
    
    def smooth_boxcar(self, width):
        """ Uses boxcar smoothing to smooth the internal arrays.
        
            Parameters
            ----------
            width : any numeric type
                The width of the smoothing in pixels (??)
        
        """
        pass
    
    def smooth_boxcar(self, width):
        """ Uses boxcar smoothing to smooth the internal arrays.
        
            Parameters
            ----------
            width : any numeric type
                The width of the smoothing in pixels (??)
        
        """
        pass
    
    def custom_smooth(self):
        """ Allows for a user defined smoothing operation.
            
            This is simply a template function meant to be overwritten
            by a custom function e.g.:
            
            def so_smooth(self):
                # Do something with flux, dispersion, error
            
            spectrumObject.custom_smooth = so_smooth
            spectrumObject.custom_smooth()
        """
        pass
    
    def plot(self, **kwargs):
        """ Plotting utility for the spectrum.
            
            Parameters
            ----------
            axes : `~matplotlib.pyplot.Axes`, optional
                A matplotlib Axes() object to plot the spectrum on
            filename : str, optional
                A filename to save the plot to
            show_error : bool, optional
                Decide whether to use the error data when plotting
                
            Notes
            -----
            Calling Spectrum1D.plot() without any arguments will attempt to
            display the object on screen using matplotlib's `show()` command.
            
            Calling Spectrum1D.plot(axes=ax) will plot the data onto the specified
            Axes() object.
            
            Calling Spectrum1D.plot(filename="myPlot.png") will plot the data onto
            a new figure and save the figure to the specified filename.
        
        """
        pass
            
    
    