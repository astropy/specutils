from __future__ import division

import numpy as np

from ..spectra import Spectrum1D, SpectralRegion

__all__ = ['excise_regions', 'linear_exciser']


def linear_exciser(spectrum, region):
    """
    Basic spectral excise method where the spectral region defined by the
    2-tuple parameter `region` (start and end wavelengths) will result
    in the flux between those regions set to a linear ramp of the
    two points immediately before and after the start and end regions.

    Other methods could be defined by the user to do other types of excision.

    Parameters
    ----------
    spectrum : `~specutils.Spectrum1D`
        The `~specutils.Spectrum1D` object to which the smoothing will be applied.

    region : `~specutils.SpectralRegion`
        Region to excise.

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
    wavelengths_in = (wavelengths >= region.lower) & (wavelengths < region.upper)
    inclusive_indices = np.nonzero(wavelengths_in)[0]

    #
    # Now set the flux values for these indices to be a
    # linear range
    #

    s, e = inclusive_indices[0]-1, inclusive_indices[-1]+1

    flux = spectrum.flux.value
    modified_flux = flux
    modified_flux[s:e] = np.linspace(flux[s], flux[e], len(inclusive_indices)+1)

    # Return a new object with the regions excised.
    return Spectrum1D(flux=modified_flux*spectrum.flux.unit,
                      spectral_axis=wavelengths,
                      uncertainty=spectrum.uncertainty,
                      wcs=spectrum.wcs, unit=spectrum.unit,
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

    regions : list of `~specutils.SpectralRegion`
        Each element of the list is a `~specutils.SpectralRegion`. The flux
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

    region : `~specutils.SpectralRegion`
        Region to excise.

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

    if not isinstance(region, SpectralRegion):
        raise ValueError('The region parameter must be a 2-tuples of start and end wavelengths.')

    #
    #  Call the exciser method
    #

    return exciser(spectrum, region)
