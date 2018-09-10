import numpy as np
from astropy.stats.funcs import gaussian_fwhm_to_sigma
from ..spectra import SpectralRegion


__all__ = ['sigma_full_width', 'full_width_half_max']


def _computation_wrapper(func, spectrum, region):

    # No region, therefore whole spectrum.
    if region is None:
        return func(spectrum)

    # Single region
    elif isinstance(region, SpectralRegion):
        return func(spectrum, region=region)

    # List of regions
    elif isinstance(region, list):
        return [func(spectrum, region=reg) for reg in region]


def sigma_full_width(spectrum, region=None):
    """
    Calculate the full width of the spectrum based on an approximate value of
    sigma.  This will be calculated over the regions, if they are specified.

    Parameters
    ----------
    spectrum : `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object overwhich the equivalent width will be calculated.

    region: `~specutils.utils.SpectralRegion` or list of `~specutils.utils.SpectralRegion`
        Region within the spectrum to calculate the gaussian sigma width.

    Returns
    -------
    full_width : `~astropy.units.Quantity` or list (based on region input)
        Approximate full width of the signal
    """
    return _computation_wrapper(_compute_sigma_full_width, spectrum, region)


def full_width_half_max(spectrum, region=None):
    """
    Calculate the full width half max of the spectrum.  This will be calculated
    over the regions, if they are specified.

    Parameters
    ----------
    spectrum : `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object overwhich the equivalent width will be calculated.

    region: `~specutils.utils.SpectralRegion` or list of `~specutils.utils.SpectralRegion`
        Region within the spectrum to calculate the FWHM value.

    Returns
    -------
    full_width_half_max : `~astropy.units.Quantity` or list (based on region input)
        Approximate full width of the signal at half max
    """
    return _computation_wrapper(_compute_full_width_half_max, spectrum, region)


def _compute_full_width_half_max(spectrum, region=None):
    """
    Calculate the full width of the spectrum at half max.

    Parameters
    ----------
    spectrum : `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object overwhich the equivalent width will be calculated.

    region: `~specutils.utils.SpectralRegion`
        Region within the spectrum to calculate the FWHM value.

    Returns
    -------
    full_width_half_max : `~astropy.units.Quantity` or list (based on region input)
        Approximate full width of the signal at half max

    Notes
    -----
    This is a helper function for the above `full_width_half_max()` method.

    """

    if region is not None:
        calc_spectrum = region.extract(spectrum)
    else:
        calc_spectrum = spectrum

    flux = calc_spectrum.flux
    spectrum = calc_spectrum.spectral_axis

    dx = spectrum - np.mean(spectrum)
    fwhm = 2 * np.sqrt(np.sum((dx * dx) * flux, axis=-1) / np.sum(flux, axis=-1))

    return fwhm


def _compute_sigma_full_width(spectrum, region=None):
    """
    Calculate the full width of the spectrum based on an approximate value of
    sigma.

    Parameters
    ----------
    spectrum : `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object overwhich the equivalent width will be calculated.

    region: `~specutils.utils.SpectralRegion`
        Region within the spectrum to calculate the gaussian sigma width.

    Returns
    -------
    full_width : `~astropy.units.Quantity` or list (based on region input)
        Approximate full width of the signal

    Notes
    -----
    This is a helper function for the above `sigma_full_width()` method.

    """

    fwhm = _compute_full_width_half_max(spectrum, region)
    sigma = fwhm * gaussian_fwhm_to_sigma

    return sigma * 2
