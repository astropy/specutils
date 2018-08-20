import numpy as np
from astropy.stats.funcs import gaussian_fwhm_to_sigma
from ..spectra import SpectralRegion

__all__ = ['sigma_full_width']


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
    full_width : float or list (based on region input)
        Approximate full width of the signal

    Notes
    -----
    The spectrum will need to have the uncertainty defined in order
    for the gaussian sigma width to be calculated.

    """

    # No region, therefore whole spectrum.
    if region is None:
        return _compute_sigma_full_width(spectrum)

    # Single region
    elif isinstance(region, SpectralRegion):
        return _compute_sigma_full_width(spectrum, region=region)

    # List of regions
    elif isinstance(region, list):
        return [_compute_sigma_full_width(spectrum, region=reg) for reg in region]


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
    full_width : float or list (based on region input)
        Approximate full width of the signal

    Notes
    -----
    This is a helper function for the above `sigma_full_width()` method.

    """

    if region is not None:
        calc_spectrum = region.extract(spectrum)
    else:
        calc_spectrum = spectrum

    flux = calc_spectrum.flux
    frequencies = calc_spectrum.frequency

    dx = frequencies - np.mean(frequencies)
    fwhm = 2 * np.sqrt(np.sum((dx * dx) * flux, axis=-1) / np.sum(flux, axis=-1))
    sigma = fwhm * gaussian_fwhm_to_sigma

    return sigma * 2
