from __future__ import division

# TODO: This can be removed when specutils has a deepcopy method.
from copy import deepcopy

from astropy import convolution
from scipy.signal import medfilt
from ..spectra import Spectrum1D

__all__ = ['convolution_smooth', 'box_smooth', 'gaussian_smooth',
           'trapezoid_smooth', 'median_smooth']


def convolution_smooth(spectrum, kernel, inplace=False):
    """
    Apply a convolution based smoothing to the spectrum. The kernel must be one
    of the 1D kernels defined in ``astropy.convolution``.

    This method can be used along but also is used by other specific methods below.

    Parameters
    ----------
    spectrum : `~specutils.spectra.Spectrum1D`
        The `~specutils.spectra.Spectrum1D` object to which the smoothing will be applied.
    kernel : ``astropy.convolution.Kernel1D`` subclass.
        The convolution based smoothing kernel.
    inplace : bool
        Choose if the smoothing should be applied to the input spectrum or a copy.
        Default is False.

    Returns
    -------
    spectrum : `~specutils.spectra.Spectrum1D`
        Output `~specutils.spectra.Spectrum1D` which is either the one passed in (inplace=True) or
        a copy of the one passed in with the updated flux (inplace=False)

    Raises
    ------
    ValueError
       In the case that ``spectrum`` and ``kernel`` are not the correct types.

    """
    # Parameter checks
    if not isinstance(spectrum, Spectrum1D):
        raise ValueError('The spectrum parameter must be a Spectrum1D object')

    if not isinstance(kernel, convolution.Kernel1D):
        raise ValueError('The kernel parameter must be a subclass of the astropy.convolution.Kernel1D')

    # Get the flux of the input spectrum
    flux = spectrum.flux

    # Smooth based on the input kernel
    smoothed_flux = convolution.convolve(flux, kernel)

    # Set the smoothed flux to the input spectrum (inplace) or
    # to a copy and return (not inplace)
#     if inplace:
#         spectrum.flux = smoothed_flux
#         return spectrum
#     else:
#         # TODO: This can be modified when specutils has a deepcopy method.
#         spectrum_out = deepcopy(spectrum)
#         spectrum_out.flux = smoothed_flux
#         return spectrum_out
    return Spectrum1D(flux=smoothed_flux, spectral_axis=spectrum.spectral_axis,
                      wcs=spectrum.wcs, unit=spectrum.unit, 
                      spectral_axis_unit=spectrum.spectral_axis_unit, 
                      velocity_convention=spectrum.velocity_convention,
                      rest_value=spectrum.rest_value)


def box_smooth(spectrum, width, inplace=False):
    """
    Smooth a `~specutils.spectra.Spectrum1D` instance based on a ``astropy.convolution.Box1D`` kernel.

    Parameters
    ----------
    spectrum : `~specutils.spectra.Spectrum1D`
        The spectrum object to which the smoothing will be applied.
    width : number
        The width of the kernel as defined in ``astropy.convolution.Box1DKernel``
    inplace : bool
        Choose if the smoothing should be applied to the input spectrum or a copy.
        Default is False.

    Returns
    -------
    spectrum : `~specutils.spectra.Spectrum1D`
        Output `~specutils.spectra.Spectrum1D` which is either the one passed in (inplace=True) or
        a copy of the one passed in with the updated flux (inplace=False)

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
    return convolution_smooth(spectrum, box1d_kernel, inplace)


def gaussian_smooth(spectrum, stddev, inplace=False):
    """
    Smooth a `~specutils.spectra.Spectrum1D` instance based on a ``astropy.convolution.Gaussian1DKernel``.

    Parameters
    ----------
    spectrum : `~specutils.spectra.Spectrum1D`
        The spectrum object to which the smoothing will be applied.
    stddev : number
        The stddev of the kernel as defined in ``astropy.convolution.Gaussian1DKernel``
    inplace : bool
        Choose if the smoothing should be applied to the input spectrum or a copy.
        Default is False.

    Returns
    -------
    spectrum : `~specutils.spectra.Spectrum1D`
        Output `~specutils.spectra.Spectrum1D` which is either the one passed in (inplace=True) or
        a copy of the one passed in with the updated flux (inplace=False)

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
    return convolution_smooth(spectrum, gaussian_kernel, inplace)


def trapezoid_smooth(spectrum, width, inplace=False):
    """
    Smoothing based on a Trapezoid kernel.

    Parameters
    ----------
    spectrum : `~specutils.spectra.Spectrum1D`
        The `~specutils.spectra.Spectrum1D` object to which the smoothing will be applied.
    width : number
        The width of the kernel as defined in ``astropy.convolution.Trapezoid1DKernel``
    inplace : bool
        Choose if the smoothing should be applied to the input spectrum or a copy.
        Default is False.

    Returns
    -------
    spectrum : `~specutils.spectra.Spectrum1D`
        Output `~specutils.spectra.Spectrum1D` which is either the one passed in (inplace=True) or
        a copy of the one passed in with the updated flux (inplace=False)

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
    return convolution_smooth(spectrum, trapezoid_kernel, inplace)


def median_smooth(spectrum, width, inplace=False):
    """
    Smoothing based on a median filter. The median filter smoothing
    is implemented using the ``scipy.signals.medfilt`` function.

    Parameters
    ----------
    spectrum : `~specutils.spectra.Spectrum1D`
        The `~specutils.spectra.Spectrum1D` object to which the smoothing will be applied.
    width : number
        The width of the median filter.
    inplace : bool
        Choose if the smoothing should be applied to the input spectrum or a copy.
        Default is False.

    Returns
    -------
    spectrum : `~specutils.spectra.Spectrum1D`
        Output `~specutils.spectra.Spectrum1D` which is either the one passed in (inplace=True) or
        a copy of the one passed in with the updated flux (inplace=False)

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

#    # Set the smoothed flux to the input spectrum (inplace) or
#    # to a copy and return (not inplace)
#    if inplace:
#        spectrum.flux = smoothed_flux
#        return spectrum
#    else:
#        # TODO: This can be modified when specutils has a deepcopy method.
#        spectrum_out = deepcopy(spectrum)
#        spectrum_out.flux = smoothed_flux
#        return spectrum_out

    return Spectrum1D(flux=smoothed_flux, spectral_axis=spectrum.spectral_axis,
                      wcs=spectrum.wcs, unit=spectrum.unit, 
                      spectral_axis_unit=spectrum.spectral_axis_unit, 
                      velocity_convention=spectrum.velocity_convention,
                      rest_value=spectrum.rest_value)
