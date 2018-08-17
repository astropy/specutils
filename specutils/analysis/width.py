import numpy as np
from astropy.stats.funcs import gaussian_fwhm_to_sigma
from ..spectra import SpectralRegion

__all__ = ['sigma']


def sigma(spectrum, region=None):
    """
    Calculate the gaussian sigma width of the spectrum.  This will be
    calculated over the regions, if they are specified.

    Parameters
    ----------
    spectrum : `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object overwhich the equivalent width will be calculated.

    region: `~specutils.utils.SpectralRegion` or list of `~specutils.utils.SpectralRegion`
        Region within the spectrum to calculate the gaussian sigma width.

    Returns
    -------
    snr : float or list (based on region input)
        Signal to noise ratio of the spectrum or within the regions

    Notes
    -----
    The spectrum will need to have the uncertainty defined in order
    for the gaussian sigma width to be calculated.

    """

    # No region, therefore whole spectrum.
    if region is None:
        return _compute_sigma(spectrum)

    # Single region
    elif isinstance(region, SpectralRegion):
        return _compute_sigma(spectrum, region=region)

    # List of regions
    elif isinstance(region, list):
        return [_compute_sigma(spectrum, region=reg) for reg in region]


def _compute_sigma(spectrum, region=None):
    """
    Calculate the gaussian sigma width of the spectrum.

    Parameters
    ----------
    spectrum : `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object overwhich the equivalent width will be calculated.

    region: `~specutils.utils.SpectralRegion`
        Region within the spectrum to calculate the gaussian sigma width.

    Returns
    -------
    snr : float or list (based on region input)
        Signal to noise ratio of the spectrum or within the regions

    Notes
    -----
    This is a helper function for the above `snr()` method.

    """

    if region is not None:
        calc_spectrum = region.extract(spectrum)
    else:
        calc_spectrum = spectrum

    flux = calc_spectrum.flux
    frequencies = calc_spectrum.frequency

    dx = frequencies - np.mean(frequencies)
    fwhm = 2 * np.sqrt(np.sum((dx * dx) * flux) / np.sum(flux))
    sigma = fwhm * gaussian_fwhm_to_sigma

    return sigma
