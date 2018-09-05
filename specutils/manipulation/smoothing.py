from __future__ import division

import astropy.units as u
from astropy import convolution
from scipy.signal import medfilt

from ..spectra import Spectrum1D

__all__ = ['convolution_smooth', 'box_smooth', 'gaussian_smooth',
           'trapezoid_smooth', 'median_smooth']


def convolution_smooth(spectrum, kernel):
    """
    Apply a convolution based smoothing to the spectrum. The kernel must be one
    of the 1D kernels defined in `astropy.convolution`.

    This method can be used along but also is used by other specific methods below.

    Parameters
    ----------
    spectrum : `~specutils.Spectrum1D`
        The `~specutils.Spectrum1D` object to which the smoothing will be applied.
    kernel : `astropy.convolution.Kernel1D` subclass or array.
        The convolution based smoothing kernel - anything that `astropy.convolution.convolve` accepts.

    Returns
    -------
    spectrum : `~specutils.Spectrum1D`
        Output `~specutils.Spectrum1D` which is copy of the one passed in with the updated flux.

    Raises
    ------
    ValueError
       In the case that ``spectrum`` and ``kernel`` are not the correct types.

    """
    # Parameter checks
    if not isinstance(spectrum, Spectrum1D):
        raise ValueError('The spectrum parameter must be a Spectrum1D object')

    # Get the flux of the input spectrum
    flux = spectrum.flux

    # Smooth based on the input kernel
    smoothed_flux = convolution.convolve(flux, kernel)

    # Return a new object with the smoothed flux.
    return Spectrum1D(flux=u.Quantity(smoothed_flux, spectrum.unit),
                      spectral_axis=u.Quantity(spectrum.spectral_axis,
                                               spectrum.spectral_axis_unit),
                      wcs=spectrum.wcs,
                      velocity_convention=spectrum.velocity_convention,
                      rest_value=spectrum.rest_value)


def box_smooth(spectrum, width):
    """
    Smooth a `~specutils.Spectrum1D` instance based on a `astropy.convolution.Box1DKernel` kernel.

    Parameters
    ----------
    spectrum : `~specutils.Spectrum1D`
        The spectrum object to which the smoothing will be applied.
    width : number
        The width of the kernel, in pixels, as defined in `astropy.convolution.Box1DKernel`

    Returns
    -------
    spectrum : `~specutils.Spectrum1D`
        Output `~specutils.Spectrum1D` which a copy of the one passed in with the updated flux.

    Raises
    ------
    ValueError
       In the case that ``width`` is not the correct type or value.

    """
    # Parameter checks
    if not isinstance(width, (int, float)) or width <= 0:
        raise ValueError('The width parameter, {}, must be a number greater than 0'.format(
                width))

    # Create the gaussian kernel
    box1d_kernel = convolution.Box1DKernel(width)

    # Call and return the convolution smoothing.
    return convolution_smooth(spectrum, box1d_kernel)


def gaussian_smooth(spectrum, stddev):
    """
    Smooth a `~specutils.Spectrum1D` instance based on a `astropy.convolution.Gaussian1DKernel`.

    Parameters
    ----------
    spectrum : `~specutils.Spectrum1D`
        The spectrum object to which the smoothing will be applied.
    stddev : number
        The stddev of the kernel, in pixels, as defined in `astropy.convolution.Gaussian1DKernel`

    Returns
    -------
    spectrum : `~specutils.Spectrum1D`
        Output `~specutils.Spectrum1D` which is copy of the one passed in with the updated flux.

    Raises
    ------
    ValueError
       In the case that ``stddev`` is not the correct type or value.

    """
    # Parameter checks
    if not isinstance(stddev, (int, float)) or stddev <= 0:
        raise ValueError('The stddev parameter, {}, must be a number greater than 0'.format(
                stddev))

    # Create the gaussian kernel
    gaussian_kernel = convolution.Gaussian1DKernel(stddev)

    # Call and return the convolution smoothing.
    return convolution_smooth(spectrum, gaussian_kernel)


def trapezoid_smooth(spectrum, width):
    """
    Smoothing based on a `astropy.convolution.Trapezoid1DKernel` kernel.

    Parameters
    ----------
    spectrum : `~specutils.Spectrum1D`
        The `~specutils.Spectrum1D` object to which the smoothing will be applied.
    width : number
        The width of the kernel, in pixels, as defined in `astropy.convolution.Trapezoid1DKernel`

    Returns
    -------
    spectrum : `~specutils.Spectrum1D`
        Output `~specutils.Spectrum1D` which is copy of the one passed in with the updated flux.

    Raises
    ------
    ValueError
       In the case that ``width`` is not the correct type or value.

    """
    # Parameter checks
    if not isinstance(width, (int, float)) or width <= 0:
        raise ValueError('The stddev parameter, {}, must be a number greater than 0'.format(
                width))

    # Create the gaussian kernel
    trapezoid_kernel = convolution.Trapezoid1DKernel(width)

    # Call and return the convolution smoothing.
    return convolution_smooth(spectrum, trapezoid_kernel)


def median_smooth(spectrum, width):
    """
    Smoothing based on a median filter. The median filter smoothing
    is implemented using the `scipy.signal.medfilt` function.

    Parameters
    ----------
    spectrum : `~specutils.Spectrum1D`
        The `~specutils.Spectrum1D` object to which the smoothing will be applied.
    width : number
        The width of the median filter in pixels.

    Returns
    -------
    spectrum : `~specutils.Spectrum1D`
        Output `~specutils.Spectrum1D` which is copy of the one passed in with the updated flux.

    Raises
    ------
    ValueError
       In the case that ``spectrum`` or ``width`` are not the correct type or value.

    """
    # Parameter checks
    if not isinstance(spectrum, Spectrum1D):
        raise ValueError('The spectrum parameter must be a Spectrum1D object')

    if not isinstance(width, (int, float)) or width <= 0:
        raise ValueError('The stddev parameter, {}, must be a number greater than 0'.format(
                width))

    # Get the flux of the input spectrum
    flux = spectrum.flux

    # Smooth based on the input kernel
    smoothed_flux = medfilt(flux, width)

    # Return a new object with the smoothed flux.
    return Spectrum1D(flux=u.Quantity(smoothed_flux, spectrum.unit),
                      spectral_axis=u.Quantity(spectrum.spectral_axis,
                                               spectrum.spectral_axis_unit),
                      wcs=spectrum.wcs,
                      velocity_convention=spectrum.velocity_convention,
                      rest_value=spectrum.rest_value)
