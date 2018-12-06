"""
A module for analysis tools focused on determining the location of
spectral features.
"""

import numpy as np
from ..spectra import SpectralRegion
from ..manipulation import extract_region


__all__ = ['centroid']


def centroid(spectrum, region):
    """
    Calculate the centroid of a region, or regions, of the spectrum.

    Parameters
    ----------
    spectrum : `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object overwhich the centroid will be calculated.

    region: `~specutils.utils.SpectralRegion` or list of `~specutils.utils.SpectralRegion`
        Region within the spectrum to calculate the centroid.

    Returns
    -------
    centroid : float or list (based on region input)
        Centroid of the spectrum or within the regions

    Notes
    -----
    The spectrum will need to be continuum subtracted before calling
    this method. See the
    `analysis documentation <https://specutils.readthedocs.io/en/latest/basic_analysis.html>`_ for more information.

    """

    # No region, therefore whole spectrum.
    if region is None:
        return _centroid_single_region(spectrum)

    # Single region
    elif isinstance(region, SpectralRegion):
        return _centroid_single_region(spectrum, region=region)

    # List of regions
    elif isinstance(region, list):
        return [_centroid_single_region(spectrum, region=reg)
                for reg in region]


def _centroid_single_region(spectrum, region=None):
    """
    Calculate the centroid of the spectrum based on the flux and uncertainty
    in the spectrum.

    Parameters
    ----------
    spectrum : `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object overwhich the centroid will be calculated.

    region: `~specutils.utils.SpectralRegion`
        Region within the spectrum to calculate the centroid.

    Returns
    -------
    centroid : float or list (based on region input)
        Centroid of the spectrum or within the regions

    Notes
    -----
    This is a helper function for the above `centroid()` method.

    """

    if region is not None:
        calc_spectrum = extract_region(spectrum, region)
    else:
        calc_spectrum = spectrum

    flux = calc_spectrum.flux
    dispersion = calc_spectrum.spectral_axis

    if len(flux.shape) > 1:
        dispersion = np.tile(dispersion, [flux.shape[0], 1])

    # the axis=-1 will enable this to run on single-dispersion, single-flux
    # and single-dispersion, multiple-flux
    return np.sum(flux * dispersion, axis=-1) / np.sum(flux, axis=-1)
