import numpy as np
from ..spectra import SpectralRegion

__all__ = ['centroid']


def centroid(spectrum, region=None):
    """
    Calculate the centroid of the spectrum based on the flux and uncertainty
    in the spectrum. This will be calculated over the regions, if they
    are specified.

    Parameters
    ----------
    spectrum : `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object overwhich the equivalent width will be calculated.

    region: `~specutils.utils.SpectralRegion` or list of `~specutils.utils.SpectralRegion`
        Region within the spectrum to calculate the centroid.

    Returns
    -------
    centroid : float or list (based on region input)
        Signal to noise ratio of the spectrum or within the regions

    Notes
    -----
    The spectrum will need to be continuum subtracted before calling
    this method.

    """

    # No region, therefore whole spectrum.
    if region is None:
        return _centroid_single_region(spectrum)

    # Single region
    elif isinstance(region, SpectralRegion):
        return _ntroid(spectrum, region=region)

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
        The spectrum object overwhich the equivalent width will be calculated.

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
        calc_spectrum = region.extract(spectrum)
    else:
        calc_spectrum = spectrum

    flux = calc_spectrum.flux
    dispersion = calc_spectrum.spectral_axis

    if len(flux.shape) > 1:
        dispersion = np.tile(dispersion[np.newaxis].T, [1, flux.shape[1]])

    # the axis=0 will enable this to run on single-dispersion, single-flux
    # and single-dispersion, multiple-flux
    return np.sum(flux * dispersion, axis=0) / np.sum(flux, axis=0)
