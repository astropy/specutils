import numpy as np
import warnings

from astropy.utils.exceptions import AstropyUserWarning
from astropy.coordinates import SpectralCoord
from ..spectra import Spectrum1D, SpectralRegion

__all__ = ['excise_regions', 'linear_exciser', 'spectrum_from_model']


def true_exciser(spectrum, region, inclusive_upper=False):
    """
    Basic spectral excise method where the array elements in the spectral
    region defined by the parameter ``region`` (a `~specutils.SpectralRegion`)
    will be deleted from all applicable elements of the Spectrum1D object:
    flux, spectral_axis, mask, and uncertainty. If multiple subregions are
    defined in ``region``, all the subregions will be excised.

    Other methods could be defined by the user to do other types of excision.

    Parameters
    ----------
    spectrum : `~specutils.Spectrum1D`
        The `~specutils.Spectrum1D` object to which the excision will be applied.

    region : `~specutils.SpectralRegion`
        The region of the spectrum to remove.

    inclusive_upper : bool, optional
        Set to True to override the default pythonic slicing on the regions (inclusive
        lower bound, exclusive upper bound) to include the specified upper bound instead.

    Returns
    -------
    spectrum : `~specutils.Spectrum1D`
        Output `~specutils.Spectrum1D` with the region excised.

    Raises
    ------
    ValueError
       In the case that ``spectrum`` and ``region`` are not the correct types.

    """

    spectral_axis = spectrum.spectral_axis
    excise_indices = None

    for subregion in region:
        # Find the indices of the spectral_axis array corresponding to the subregion
        if inclusive_upper:
            region_mask = (spectral_axis >= subregion.lower) & (spectral_axis <= subregion.upper)
        else:
            region_mask = (spectral_axis >= subregion.lower) & (spectral_axis < subregion.upper)

        temp_indices = np.nonzero(region_mask)[0]
        if excise_indices is None:
            excise_indices = temp_indices
        else:
            excise_indices = np.hstack(excise_indices, temp_indices)

    new_flux = np.delete(spectrum.flux, excise_indices)
    new_spectral_axis = np.delete(spectrum.spectral_axis, excise_indices)

    if spectrum.mask is not None:
        new_mask = np.delete(spectrum.mask, excise_indices)
    else:
        new_mask = None

    if spectrum.uncertainty is not None:
        new_uncertainty = spectrum.uncertainty.__class__(
            np.delete(spectrum.uncertainty.array, excise_indices),
            unit=spectrum.uncertainty.unit)
    else:
        new_uncertainty = None

    # Return a new object with the regions excised.
    return Spectrum1D(flux=new_flux,
                      spectral_axis=new_spectral_axis,
                      uncertainty=new_uncertainty,
                      mask=new_mask,
                      wcs=spectrum.wcs,
                      velocity_convention=spectrum.velocity_convention,
                      rest_value=spectrum.rest_value if not isinstance(new_spectral_axis, SpectralCoord) else None,
                      radial_velocity=spectrum.radial_velocity if not isinstance(new_spectral_axis, SpectralCoord) else None)


def linear_exciser(spectrum, region, inclusive_upper=False):
    """
    Basic spectral excise method where the spectral region defined by the
    parameter ``region`` (a `~specutils.SpectralRegion`) will result
    in the flux between those regions set to a linear ramp of the
    two points immediately before and after the start and end of the region.

    Other methods could be defined by the user to do other types of excision.

    Parameters
    ----------
    spectrum : `~specutils.Spectrum1D`
        The `~specutils.Spectrum1D` object to which the excision will be applied.

    region : `~specutils.SpectralRegion`
        The region of the spectrum to replace.

    Returns
    -------
    spectrum : `~specutils.Spectrum1D`
        Output `~specutils.Spectrum1D` with the region excised.

    Raises
    ------
    ValueError
       In the case that ``spectrum`` and ``region`` are not the correct types.

    """

    spectral_axis = spectrum.spectral_axis.copy()
    flux = spectrum.flux.copy()
    modified_flux = flux
    if spectrum.uncertainty is not None:
        new_uncertainty = spectrum.uncertainty.copy()
    else:
        new_uncertainty = None

    # Need to add a check that the subregions don't overlap, since that could
    # cause undesired results. For now warn if there is more than one subregion
    if len(region) > 1:
        # Raise a warning if the SpectralRegion has more than one subregion, since
        # the handling for this is perhaps unexpected
        warnings.warn("A SpectralRegion with multiple subregions was provided as "
             "input. This may lead to undesired behavior with linear_exciser if "
             "the subregions overlap.",
             AstropyUserWarning)

    for subregion in region:
        # Find the indices of the spectral_axis array corresponding to the subregion
        if inclusive_upper:
            region_mask = (spectral_axis >= subregion.lower) & (spectral_axis <= subregion.upper)
        else:
            region_mask = (spectral_axis >= subregion.lower) & (spectral_axis < subregion.upper)

        inclusive_indices = np.nonzero(region_mask)[0]
        # Now set the flux values for these indices to be a
        # linear range
        s, e = max(inclusive_indices[0]-1, 0), min(inclusive_indices[-1]+1,
                                               spectral_axis.size-1)

        modified_flux[s:e+1] = np.linspace(flux[s], flux[e], modified_flux[s:e+1].size)

        # Add the uncertainty of the two linear interpolation endpoints in
        # quadrature and apply to the excised region.
        if new_uncertainty is not None:
            new_uncertainty[s:e] = np.sqrt(spectrum.uncertainty[s]**2 + spectrum.uncertainty[e]**2)

    # Return a new object with the regions excised.
    return Spectrum1D(flux=modified_flux,
                      spectral_axis=spectral_axis,
                      uncertainty=new_uncertainty,
                      wcs=spectrum.wcs,
                      mask=spectrum.mask,
                      velocity_convention=spectrum.velocity_convention,
                      rest_value=spectrum.rest_value if not isinstance(spectral_axis, SpectralCoord) else None,
                      radial_velocity=spectrum.radial_velocity if not isinstance(spectral_axis, SpectralCoord) else None)


def excise_regions(spectrum, regions, exciser=true_exciser, inclusive_upper=False):
    """
    Method to remove or replace the flux in the defined regions of the spectrum
    depending on the function provided in the ``exciser`` argument.

    Parameters
    ----------
    spectrum : `~specutils.Spectrum1D`
        The `~specutils.Spectrum1D` object to which the excision will be applied.

    regions : list of `~specutils.SpectralRegion`
        Each element of the list is a `~specutils.SpectralRegion`. The flux
        between the lower and upper spectral axis value of each region will be "cut out"
        and optionally replaced with interpolated values using the ``exciser`` method.
        Note that non-overlapping regions should be provided as separate
        `~specutils.SpectralRegion` objects in this list, not as sub-regions
        in a single object in the list.

    exciser : function
        Method that takes the spectrum and region and does the excising. Other
        methods could be defined and used by this routine.
        default: true_exciser

    inclusive_upper : bool, optional
        Set to True to override the default pythonic slicing on the regions (inclusive
        lower bound, exclusive upper bound) to include the specified upper bound instead.

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
        spectrum = excise_region(spectrum, region, exciser, inclusive_upper=inclusive_upper)

    return spectrum


def excise_region(spectrum, region, exciser=true_exciser, inclusive_upper=False):
    """
    Method to remove or replace the flux in the defined regions of the spectrum
    depending on the function provided in the ``exciser`` argument.

    Parameters
    ----------
    spectrum : `~specutils.Spectrum1D`
        The `~specutils.Spectrum1D` object to which the smoothing will be applied.

    region : `~specutils.SpectralRegion`
        A `~specutils.SpectralRegion` object defining the region to excise.
        If excising multiple regions is desired, they should be input as a
        list of separate `~specutils.SpectralRegion` objects to
        ``excise_regions``, not as subregions defined in a single
        `~specutils.SpectralRegion`.

    exciser: method
        Method that takes the spectrum and region and does the excising. Other
        methods could be defined and used by this routine.
        default: true_exciser

    inclusive_upper : bool, optional
        Set to True to override the default pythonic slicing on the regions (inclusive
        lower bound, exclusive upper bound) to include the specified upper bound instead.

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
        raise ValueError('The spectrum parameter must be a Spectrum1D object.')

    if not isinstance(region, SpectralRegion):
        raise ValueError('The region parameter must be a SpectralRegion object.')

    #
    #  Call the exciser method
    #

    return exciser(spectrum, region, inclusive_upper=inclusive_upper)


def spectrum_from_model(model_input, spectrum):
    """
    This method will create a `~specutils.Spectrum1D` object
    with the flux defined by calling the input ``model``. All
    other parameters for the output `~specutils.Spectrum1D` object
    will be the same as the input `~specutils.Spectrum1D` object.

    Parameters
    ----------
    model : `~astropy.modeling.Model`
        The input model or compound model from which flux is calculated.

    spectrum : `~specutils.Spectrum1D`
        The `~specutils.Spectrum1D` object to use as the model template.

    Returns
    -------
    spectrum : `~specutils.Spectrum1D`
        Output `~specutils.Spectrum1D` which is copy of the one passed in with the updated flux.
        The uncertainty will not be copied as it is not necessarily the same.

    """

    # If the input model has units then we will call it normally.
    if getattr(model_input, model_input.param_names[0]).unit is not None:
        flux = model_input(spectrum.spectral_axis)

    # If the input model does not have units, then assume it is in
    # the same units as the input spectrum.
    else:
        flux = model_input(spectrum.spectral_axis.value)*spectrum.flux.unit

    return Spectrum1D(flux=flux,
                      spectral_axis=spectrum.spectral_axis,
                      wcs=spectrum.wcs,
                      velocity_convention=spectrum.velocity_convention,
                      rest_value=spectrum.rest_value)
