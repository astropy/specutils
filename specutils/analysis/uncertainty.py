"""
A module for analysis tools dealing with uncertainties or error analysis in
spectra.
"""

import numpy as np
from ..spectra import SpectralRegion
from ..manipulation import extract_region

__all__ = ['snr']


def snr(spectrum, region=None):
    """
    Calculate the mean S/N of the spectrum based on the flux and uncertainty
    in the spectrum. This will be calculated over the regions, if they
    are specified.

    Parameters
    ----------
    spectrum : `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object overwhich the equivalent width will be calculated.

    region: `~specutils.utils.SpectralRegion` or list of `~specutils.utils.SpectralRegion`
        Region within the spectrum to calculate the SNR.

    Returns
    -------
    snr : `~astropy.units.Quantity` or list (based on region input)
        Signal to noise ratio of the spectrum or within the regions

    Notes
    -----
    The spectrum will need to have the uncertainty defined in order for the SNR
    to be calculated. If the goal is instead signal to noise *per pixel*, this
    should be computed directly as ``spectrum.flux / spectrum.uncertainty``.

    """

    if not hasattr(spectrum, 'uncertainty') or spectrum.uncertainty is None:
        raise Exception("Spectrum1D currently requires the uncertainty be defined.")

    # No region, therefore whole spectrum.
    if region is None:
        return _snr_single_region(spectrum)

    # Single region
    elif isinstance(region, SpectralRegion):
        return _snr_single_region(spectrum, region=region)

    # List of regions
    elif isinstance(region, list):
        return [_snr_single_region(spectrum, region=reg)
                for reg in region]


def _snr_single_region(spectrum, region=None):
    """
    Calculate the mean S/N of the spectrum based on the flux and uncertainty
    in the spectrum.

    Parameters
    ----------
    spectrum : `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object overwhich the equivalent width will be calculated.

    region: `~specutils.utils.SpectralRegion`
        Region within the spectrum to calculate the SNR.

    Returns
    -------
    snr : `~astropy.units.Quantity` or list (based on region input)
        Signal to noise ratio of the spectrum or within the regions

    Notes
    -----
    This is a helper function for the above `snr()` method.

    """

    if region is not None:
        calc_spectrum = extract_region(spectrum, region)
    else:
        calc_spectrum = spectrum

    flux = calc_spectrum.flux
    uncertainty = calc_spectrum.uncertainty.array * spectrum.uncertainty.unit

    # the axis=-1 will enable this to run on single-dispersion, single-flux
    # and single-dispersion, multiple-flux
    return np.mean(flux / uncertainty, axis=-1)
