"""
A module for analysis tools focused on determining the moment of
spectral features.
"""

import numpy as np
from ..manipulation import extract_region
from .utils import computation_wrapper


__all__ = ['moment']


def moment(spectrum, regions=None, order=0, axis=-1):
    """
    Estimate the moment of the spectrum.


    Parameters
    ----------
    spectrum : `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object over which the width will be calculated.

    regions: `~specutils.utils.SpectralRegion` or list of `~specutils.utils.SpectralRegion`
        Region within the spectrum to calculate the gaussian sigma width. If
        regions is `None`, computation is performed over entire spectrum.

    order : int
        The order of the moment to be calculated. Default=0

    axis : int
        Axis along which a moment is calculated. Default=-1, computes along
        the last axis (spectral axis).


    Returns
    -------
    moment: `float` or list (based on region input)
        Moment of the spectrum. Returns None if (order < 0 or None)

    """
    return computation_wrapper(_compute_moment, spectrum, regions,
                               order=order, axis=axis)


def _compute_moment(spectrum, regions=None, order=0, axis=-1):
    """
    This is a helper function for the above `moment()` method.
    """
    if regions is not None:
        calc_spectrum = extract_region(spectrum, regions)
    else:
        calc_spectrum = spectrum

    # Ignore masks for now. This should be fully addressed when
    # specutils gets revamped to handle multi-dimensional masks.
    flux = calc_spectrum.flux
    spectral_axis = calc_spectrum.spectral_axis

    if order is None or order < 0:
        return None

    if order == 0:
        return np.sum(flux, axis=axis)

    dispersion = spectral_axis
    if len(flux.shape) > len(spectral_axis.shape):
        _shape = flux.shape[:-1] + (1,)
        dispersion = np.tile(spectral_axis, _shape)

    if order == 1:
        return np.sum(flux * dispersion, axis=axis) / np.sum(flux, axis=axis)

    if order > 1:
        m0 = np.sum(flux, axis=axis)
        m1 = np.sum(flux * dispersion, axis=axis) / np.sum(flux, axis=axis)

        if len(flux.shape) > 1 and (axis == len(flux.shape)-1 or axis == -1):
            _shape = flux.shape[-1:] + tuple(np.ones(flux.ndim - 1, dtype='i'))
            m1 = np.tile(m1, _shape).T

        return np.sum(flux * (dispersion - m1) ** order, axis=axis) / m0
