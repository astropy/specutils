from __future__ import division

from astropy import convolution


__all__ = ['convolution_smooth', 'box_smooth', 'gaussian_smooth',
           'mexicanhat_smooth', 'trapezoid_smooth']


def convolution_smooth(spectrum, kernel, inplace=False):
    """
    Does a naive equivalent width measures on the spectrum object.

    Parameters
    ----------
    spectrum : Spectrum1D
        The spectrum object overwhich the equivalent width will be calculated.
    kernel : astropy.modeling.models.Fittable1DModel
        The convolution based smoothing kernel.
    inplace : bool
        Choose if the smoothing should be applied to the input spectrum or a copy.

    Returns
    -------
    spectrum : spectrum1D
        Equivalent width calculation.

    TODO:  what frame of reference do you want the spectral_axis to be in ???
    """

    # Get the flux of the input spectrum
    flux = spectrum.flux

    # Smooth based on the input kernel
    smoothed_flux = convolution.convolve(flux, kernel)

    # Set the smoothed flux to the input spectrum (inplace) or
    # to a copy and return (not inplace)
    if inplace:
        spectrum.flux = smoothed_flux
        return spectrum
    else:
        spectrum_out = spectrum.copy()
        spectrum_out.flux = smoothed_flux
        return spectrum_out


def box_smooth(spectrum, width, inplace=False):
    """
    Smoothing based on a Gaussian kernel.

    Parameters
    ----------
    spectrum : Spectrum1D
        The spectrum object overwhich the equivalent width will be calculated.
    width : number
        The width of the kernel as defined in astropy.convolution.Box1DKernel
    inplace : bool
        Choose if the smoothing should be applied to the input spectrum or a copy.

    Returns
    -------
    spectrum : spectrum1D
        Equivalent width calculation.
    """
    # Create the gaussian kernel
    box1d_kernel = convolution.Box1D(width)

    # Call and return the convolution smoothing.
    return convolution_smooth(spectrum, box1d_kernel, inplace)


def gaussian_smooth(spectrum, stddev, inplace=False):
    """
    Smoothing based on a Gaussian kernel.

    Parameters
    ----------
    spectrum : Spectrum1D
        The spectrum object overwhich the equivalent width will be calculated.
    stddev : number
        The stddev of the kernel as defined in astropy.convolution.Gaussian1DKernel
    inplace : bool
        Choose if the smoothing should be applied to the input spectrum or a copy.

    Returns
    -------
    spectrum : spectrum1D
        Equivalent width calculation.
    """
    # Create the gaussian kernel
    gaussian_kernel = convolution.Gaussian1D(stddev)

    # Call and return the convolution smoothing.
    return convolution_smooth(spectrum, gaussian_kernel, inplace)


def mexicanhat_smooth(spectrum, width, inplace=False):
    """
    Smoothing based on a Gaussian kernel.

    Parameters
    ----------
    spectrum : Spectrum1D
        The spectrum object overwhich the equivalent width will be calculated.
    width : number
        The width of the kernel as defined in astropy.convolution.MexicanHat1DKernel
    inplace : bool
        Choose if the smoothing should be applied to the input spectrum or a copy.

    Returns
    -------
    spectrum : spectrum1D
        Equivalent width calculation.
    """
    # Create the gaussian kernel
    mexicanhat_kernel = convolution.MexicanHat1D(width)

    # Call and return the convolution smoothing.
    return convolution_smooth(spectrum, mexicanhat_kernel, inplace)


def trapezoid_smooth(spectrum, width, inplace=False):
    """
    Smoothing based on a Gaussian kernel.

    Parameters
    ----------
    spectrum : Spectrum1D
        The spectrum object overwhich the equivalent width will be calculated.
    width : number
        The width of the kernel as defined in astropy.convolution.Trapezoid1DKernel
    inplace : bool
        Choose if the smoothing should be applied to the input spectrum or a copy.

    Returns
    -------
    spectrum : spectrum1D
        Equivalent width calculation.
    """
    # Create the gaussian kernel
    trapezoid_kernel = convolution.Trapezoid1D(width)

    # Call and return the convolution smoothing.
    return convolution_smooth(spectrum, trapezoid_kernel, inplace)
