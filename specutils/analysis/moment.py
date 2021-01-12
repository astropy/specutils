"""
A module for analysis tools focused on determining the moment of
spectral features.
"""

import numpy as np
from ..manipulation import extract_region
from .utils import computation_wrapper


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
        Moment of the spectrum. Returns None if (order < 0 or None)

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

    if order is None or order < 0:
        return None

    if order == 0:
        # the axis=-1 will enable this to run on single-dispersion, single-flux
        # and single-dispersion, multiple-flux
        return np.sum(flux, axis=-1)

    dispersion = spectral_axis
    if len(flux.shape) > 1:
        dispersion = np.tile(spectral_axis, [flux.shape[0], 1])

    if order == 1:
        return np.sum(flux * dispersion, axis=-1) / np.sum(flux, axis=-1)

    if order > 1:
        m0 = np.sum(flux, axis=-1)
        m1 = np.sum(flux * dispersion, axis=-1) / np.sum(flux, axis=-1)

        return np.sum(flux * (dispersion - m1)**order, axis=-1) / m0
