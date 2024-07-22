import copy
import warnings

import astropy.units as u
import numpy as np
from astropy import convolution
from astropy.nddata import StdDevUncertainty, VarianceUncertainty, InverseVariance
from astropy.utils.exceptions import AstropyUserWarning
from scipy.signal import medfilt

from ..spectra import Spectrum1D

__all__ = ['convolution_smooth', 'box_smooth', 'gaussian_smooth',
           'trapezoid_smooth', 'median_smooth']


def convolution_smooth(spectrum, kernel):
    """
    Apply a convolution based smoothing to the spectrum. The kernel must be one
    of the 1D kernels defined in `astropy.convolution`, and will be applied along
    the spectral axis of the flux.

    This method can be used alone but also is used by other specific methods
    below.

    If the spectrum uncertainty exists and is ``StdDevUncertainty``,
    ``VarianceUncertainty`` or ``InverseVariance`` then the errors will be
    propagated through the convolution using a standard propagation of errors.
    The covariance is not considered, currently.

    Parameters
    ----------
    spectrum : `~specutils.Spectrum1D`
        The `~specutils.Spectrum1D` object to which the smoothing will be
        applied.
    kernel : `astropy.convolution.Kernel1D` subclass or array.
        The convolution based smoothing kernel - anything that
        `astropy.convolution.convolve` accepts.

    Returns
    -------
    spectrum : `~specutils.Spectrum1D`
        Output `~specutils.Spectrum1D` which is copy of the one passed in with
        the updated flux.

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

    # Expand kernel with empty leading dimensions if flux is multidimensional
    # and kernel is 1D.
    if isinstance(kernel, np.ndarray):
        kernel_ndim = kernel.ndim
    else:
        kernel_ndim = kernel.array.ndim

    if flux.ndim > 1 and kernel_ndim == 1:
        expand_axes = tuple(np.arange(flux.ndim-1))
        kernel = np.expand_dims(kernel, expand_axes)

    # Smooth based on the input kernel
    smoothed_flux = convolution.convolve(flux, kernel)

    # Propagate the uncertainty if it exists...
    uncertainty = copy.deepcopy(spectrum.uncertainty)
    if uncertainty is not None:
        if isinstance(uncertainty, StdDevUncertainty):
            # Convert
            values = uncertainty.array
            ivar_values = 1 / values**2

            # Propagate
            prop_ivar_values = convolution.convolve(ivar_values, kernel)

            # Put back in
            uncertainty.array = 1 / np.sqrt(prop_ivar_values)

        elif isinstance(uncertainty, VarianceUncertainty):
            # Convert
            values = uncertainty.array
            ivar_values = 1 / values

            # Propagate
            prop_ivar_values = convolution.convolve(ivar_values, kernel)

            # Put back in
            uncertainty.array = 1 / prop_ivar_values

        elif isinstance(uncertainty, InverseVariance):
            # Convert
            ivar_values = uncertainty.array

            # Propagate
            prop_ivar_values = convolution.convolve(ivar_values, kernel)

            # Put back in
            uncertainty.array = prop_ivar_values
        else:
            uncertainty = None
            warnings.warn(
                "Uncertainty is {} but convolutional error propagation is "
                "not defined for that type. Uncertainty will be dropped in "
                "the convolved spectrum.".format(type(uncertainty)),
                AstropyUserWarning)

    # Return a new object with the smoothed flux.
    return spectrum._copy(flux=u.Quantity(smoothed_flux, spectrum.unit),
                          spectral_axis=u.Quantity(spectrum.spectral_axis,
                                                   spectrum.spectral_axis.unit),
                          uncertainty=uncertainty)


def box_smooth(spectrum, width):
    """
    Smooth a `~specutils.Spectrum1D` instance along the spectral axis
    based on a `astropy.convolution.Box1DKernel` kernel.

    Parameters
    ----------
    spectrum : `~specutils.Spectrum1D`
        The spectrum object to which the smoothing will be applied.
    width : number
        The width of the kernel, in pixels, as defined in
        `astropy.convolution.Box1DKernel`

    Returns
    -------
    spectrum : `~specutils.Spectrum1D`
        Output `~specutils.Spectrum1D` which a copy of the one passed in with
        the updated flux.

    Raises
    ------
    ValueError
       In the case that ``width`` is not the correct type or value.

    """
    # Parameter checks
    if not isinstance(width, (int, float)) or width <= 0:
        raise ValueError("The width parameter, {}, must be a number greater "
                         "than 0".format(width))

    # Create the gaussian kernel
    box1d_kernel = convolution.Box1DKernel(width)

    # Call and return the convolution smoothing.
    return convolution_smooth(spectrum, box1d_kernel)


def gaussian_smooth(spectrum, stddev):
    """
    Smooth a `~specutils.Spectrum1D` instance along the spectral axis
    based on a `astropy.convolution.Gaussian1DKernel`.

    Parameters
    ----------
    spectrum : `~specutils.Spectrum1D`
        The spectrum object to which the smoothing will be applied.
    stddev : number
        The stddev of the kernel, in pixels, as defined in
        `astropy.convolution.Gaussian1DKernel`

    Returns
    -------
    spectrum : `~specutils.Spectrum1D`
        Output `~specutils.Spectrum1D` which is copy of the one passed in with
        the updated flux.

    Raises
    ------
    ValueError
       In the case that ``stddev`` is not the correct type or value.

    """
    # Parameter checks
    if not isinstance(stddev, (int, float)) or stddev <= 0:
        raise ValueError("The stddev parameter, {}, must be a number greater "
                         "than 0".format(stddev))

    # Create the gaussian kernel
    gaussian_kernel = convolution.Gaussian1DKernel(stddev)

    # Call and return the convolution smoothing.
    return convolution_smooth(spectrum, gaussian_kernel)


def trapezoid_smooth(spectrum, width):
    """
    Smooth a `~specutils.Spectrum1D` instance along the spectral axis
    based on a `astropy.convolution.Trapezoid1DKernel` kernel.

    Parameters
    ----------
    spectrum : `~specutils.Spectrum1D`
        The `~specutils.Spectrum1D` object to which the smoothing will be
        applied.
    width : number
        The width of the kernel, in pixels, as defined in
        `astropy.convolution.Trapezoid1DKernel`

    Returns
    -------
    spectrum : `~specutils.Spectrum1D`
        Output `~specutils.Spectrum1D` which is copy of the one passed in with
        the updated flux.

    Raises
    ------
    ValueError
       In the case that ``width`` is not the correct type or value.

    """
    # Parameter checks
    if not isinstance(width, (int, float)) or width <= 0:
        raise ValueError("The stddev parameter, {}, must be a number greater "
                         "than 0".format(width))

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
        The `~specutils.Spectrum1D` object to which the smoothing will be
        applied.
    width : number
        The width of the median filter in pixels.

    Returns
    -------
    spectrum : `~specutils.Spectrum1D`
        Output `~specutils.Spectrum1D` which is copy of the one passed in with
        the updated flux.

    Raises
    ------
    ValueError
       In the case that ``spectrum`` or ``width`` are not the correct type or
       value.

    """
    # Parameter checks
    if not isinstance(spectrum, Spectrum1D):
        raise ValueError('The spectrum parameter must be a Spectrum1D object')

    if not isinstance(width, (int, float)) or width <= 0:
        raise ValueError("The stddev parameter, {}, must be a number greater "
                         "than 0".format(width))

    # Get the flux of the input spectrum
    flux = spectrum.flux

    if (
            not isinstance(flux.dtype, (float, int)) or
            not np.issubdtype(flux.dtype, (np.floating, np.integer))
    ):
        flux = flux.astype(float)

    # Smooth based on the input kernel
    smoothed_flux = medfilt(flux, width)

    # Return a new object with the smoothed flux.
    return spectrum._copy(flux=u.Quantity(smoothed_flux, spectrum.unit),
                          spectral_axis=u.Quantity(spectrum.spectral_axis,
                                                   spectrum.spectral_axis.unit))
