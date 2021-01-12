"""
A module for analysis tools focused on determining the moment of
spectral features.
"""

import numpy as np
from astropy.stats.funcs import gaussian_sigma_to_fwhm
from astropy.modeling.models import Gaussian1D
from ..manipulation import extract_region
from . import centroid
from .utils import computation_wrapper
from scipy.signal import chirp, find_peaks, peak_widths


__all__ = ['moment']


def moment(spectrum, regions=None, order=0):
    """
    Estimate the moment of the spectrum.


    Parameters
    ----------
    spectrum : `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object over which the width will be calculated.

    regions: `~specutils.utils.SpectralRegion` or list of `~specutils.utils.SpectralRegion`
        Region within the spectrum to calculate the gaussian sigma width. If
        regions is `None`, computation is performed over entire spectrum.

    order : int, default: 0
        The order of the moment to be calculated.

    Returns
    -------
    moment: `float` or list (based on region input)
        Moment of the spectrum

    """
    return computation_wrapper(_compute_moment, spectrum, regions, order=order)


def _compute_moment(spectrum, regions=None, order=0):
    """
    This is a helper function for the above `gaussian_sigma_width()` method.
    """

    if regions is not None:
        calc_spectrum = extract_region(spectrum, regions)
    else:
        calc_spectrum = spectrum

    if hasattr(spectrum, 'mask') and spectrum.mask is not None:
        flux = calc_spectrum.flux[~spectrum.mask]
        spectral_axis = calc_spectrum.spectral_axis[~spectrum.mask]
    else:
        flux = calc_spectrum.flux
        spectral_axis = calc_spectrum.spectral_axis

    centroid_result = centroid(spectrum, regions)

    if flux.ndim > 1:
        spectral_axis = np.broadcast_to(spectral_axis, flux.shape, subok=True)
        centroid_result = centroid_result[:, np.newaxis]

    dx = (spectral_axis - centroid_result)
    sigma = np.sqrt(np.sum((dx * dx) * flux, axis=-1) / np.sum(flux, axis=-1))

    return sigma
