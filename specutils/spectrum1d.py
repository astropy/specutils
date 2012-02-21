# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module implements the Spectrum1D class. It is eventually supposed to migrate to astropy core

# Python packages
from abc import ABCMeta, abstractmethod, abstractproperty
from copy import copy, deepcopy

# External packages
from scipy import interpolate
import numpy as np  

# Local packages
from astropy.nddata import NDData

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
    This class is abstract and can only be subclassed and used by
    classes that will access its API.

    """
    
    # This is an abstract class and can only be subclassed.
    __metaclass__ = ABCMeta
    
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
    def __init__(self):
        pass # not yet implemented
        
    @property
    def spectrum1DBase(self):
        """ This method returns an object equivalent to this spectrum but as
            a Spectrum1DBase object, i.e. without a dispersion array.
            It is left to subclasses to handle this in more detail (e.g., see Spectrum1D).
        """
        return self    
    
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
    # This is an abstract class and can only be subclassed.
    __metaclass__ = ABCMeta
    
    @property
    def spectrum1DBase(self):
        """ Return a new Spectrum1DBase object from this spectrum.
        
            This basically is the same thing but without the dispersion array,
            but we'll see how this develops. """
        # There will be a better way to do this, but this is the general idea.
        # Basically I want to create a Spectrum1DBase object from a Spectrum1D object.
        spectrum_copy = deepcopy(self)
        spectrum_copy.dispersion = None
        return spectrum_copy

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
        
    def interpolate_bspline(self, new_dispersion):
        """ Uses B-spline interpolation to resample the internal arrays onto
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
            width : integer
                The boxcar width in pixels.
                
            Reference: <a href="http://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.boxcar.html" target="_blank">scipy.signal.boxcar</a>
        
        """
        # I don't know if this is really the proper way to do this... should be tested!!
        from scipy.signal import convolve, boxcar
        self.flux = convolve(self.flux, boxcar(M=width))
    
    def smooth_custom(self):
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
                
            Usage
            -----
            s = Spectrum1D(...)
            
            s.plot() - Attempt to display the object on screen
                       using matplotlib's `show()` command.

            s.plot(axes=ax) - Plot the data onto the specified
                              Axes() object.
            
            s.plot(axes=ax) - Plot the data onto the specified
                              Axes() object.
            
            s.plot(filename="myPlot.png") - Plot the data onto a new figure and
                                save to the specified filename and path
                                (default path if not specified).
                                The format of the file will be deduced from the
                                extension (e.g. ps, png, pdf).
            
            Notes
            -----
            Where it makes sense (i.e. not conflicting), multiple parameters
            can be specified. For example,
            
            s.plot(filename="myPlot.pdf", axes=ax)
            
            will both write a file to the disk and return a matplotlib.Axes() object.
        
        """
        pass
            
class Spectrum1DCollection(object):
    """ A collection object for spectra that share the same dispersion information.
    
    """
    def __init__(self):
        self.spectra = list()
        self.dispersion = None
        
    def append(self, spectrum):
        """ Add a spectrum (of type Spectrum1D or Spectrum1DBase) to this collection.
        
        """
        self.spectra.append(new_spectrum.spectrum1DBase)
    
    def len(self):
        return len(self.spectra)
    
    @property
    def dispersion(self):
        return self.dispersion
    
    @dispersion.setter
    def dispersion(self, new_dispersion):
        """ Set the dispersion array to be used for all spectra in this collection.
        
        The dispersion argument accepts both a numpy array and a Spectrum1D object.
        When the latter is specified, the dispersion is extracted from the object.
        """
        if is_instance(new_dispersion, Spectrum1D):
            self.dispersion = new_dispersion.dispersion
        elif is_instance(new_dispersion, list):
            self.dispersion = np.array(new_dispersion)
        elif is_instance(new_dispersion, np.array):
            self.dispersion = new_dispersion
        else:
            raise ValueError, "The dispersion specified could is not a known type. Should be a list, a numpy array, or a Spectrum1D object."
        # Add other types that unambiguously define a dispersion array.

