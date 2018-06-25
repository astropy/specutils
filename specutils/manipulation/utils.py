from __future__ import division

import numpy as np

from ..spectra import Spectrum1D

__all__ = ['excise_regions']


def linear_exciser(spectrum, region):
    """
    Basic spectral excise method where the spectral region defined by the
    2-tuple parameter `region` (start and end wavelengths) will result
    in the flux between those regions be set to a linear ramp of the
    two points immediately before and after the start and end regions.

    Other methods could be defined by the user to do other types of excision.

    Parameters
    ----------
    spectrum : `~specutils.Spectrum1D`
        The `~specutils.Spectrum1D` object to which the smoothing will be applied.

    region : 2-tuple
        Region to excise, defined as a (start wavelength, end wavelength) with
        astropy units.

    Returns
    -------
    spectrum : `~specutils.Spectrum1D`
        Output `~specutils.Spectrum1D` with the region excised.

    Raises
    ------
    ValueError
       In the case that ``spectrum`` and ``regions`` are not the correct types.

    """

    #
    # Find the indices of the wavelengths in the range `range`
    #

    wavelengths = spectrum.spectral_axis
    wavelengths_in = (wavelengths > range[0]) & (wavelengths < range[1])
    inclusive_indices = np.nonzero(wavelengths_in)[0]

    #
    # Now set the flux values for these indices to be a
    # linear range
    #

    s, e = inclusive_indices[0]-1, inclusive_indices[1]+1

    flux = spectrum.flux
    modified_flux = flux
    modified_flux[s:e] = np.linspace(flux[s], flux[e], len(inclusive_indices))

    # Return a new object with the regions excised.
    return Spectrum1D(flux=modified_flux, spectral_axis=wavelengths,
                      wcs=spectrum.wcs, unit=spectrum.unit,
                      spectral_axis_unit=spectrum.spectral_axis_unit,
                      velocity_convention=spectrum.velocity_convention,
                      rest_value=spectrum.rest_value)


def excise_regions(spectrum, regions, exciser=linear_exciser):
    """
    Apply a convolution based smoothing to the spectrum. The kernel must be one
    of the 1D kernels defined in `astropy.convolution`.

    This method can be used along but also is used by other specific methods below.

    Parameters
    ----------
    spectrum : `~specutils.Spectrum1D`
        The `~specutils.Spectrum1D` object to which the smoothing will be applied.

    regions : list of 2-tuples
        Each element of the list is a 2-tuple of wavelengths. The flux
        between these wavelengths will be "cut out" using the `exciser`
        method.

    exciser: method
        Method that takes the spectrum and region and does the excising. Other
        methods could be defined and used by this routine.
        default: linear_exciser

    Returns
    -------
    spectrum : `~specutils.Spectrum1D`
        Output `~specutils.Spectrum1D` which has the regions excised.

    Raises
    ------
    ValueError
       In the case that ``spectrum`` and ``regions`` are not the correct types.

    """

    # Parameter checks
    if not isinstance(spectrum, Spectrum1D):
        raise ValueError('The spectrum parameter must be Spectrum1D object.')

    if not isinstance(regions, list):
        raise ValueError('The regions parameter must be a list of 2-tuples.')

    for region in regions:
        spectrum = excise_region(spectrum, region, exciser)

    return spectrum


def excise_region(spectrum, region, exciser=linear_exciser):
    """
    Apply a convolution based smoothing to the spectrum. The kernel must be one
    of the 1D kernels defined in `astropy.convolution`.

    This method can be used along but also is used by other specific methods below.

    Parameters
    ----------
    spectrum : `~specutils.Spectrum1D`
        The `~specutils.Spectrum1D` object to which the smoothing will be applied.

    region : 2-tuple
        Region to excise, defined as a (start wavelength, end wavelength) in
        some astropy wavelength units.

    exciser: method
        Method that takes the spectrum and region and does the excising. Other
        methods could be defined and used by this routine.
        default: linear_exciser

    Returns
    -------
    spectrum : `~specutils.Spectrum1D`
        Output `~specutils.Spectrum1D` with the region excised.

    Raises
    ------
    ValueError
       In the case that ``spectrum`` and ``region`` are not the correct types.

    """

    # Parameter checks
    if not isinstance(spectrum, Spectrum1D):
        raise ValueError('The spectrum parameter must be Spectrum1D object.')

    if not isinstance(region, tuple) or not len(region) == 2:
        raise ValueError('The region parameter must be a 2-tuples of start and end wavelengths.')

    #
    #  Call the exciser method
    #

    return exciser(spectrum, region)
