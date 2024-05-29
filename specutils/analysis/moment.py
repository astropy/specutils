"""
A module for analysis tools focused on determining the moment of
spectral features.
"""

import numpy as np
from ..manipulation import extract_region
from ..spectra import SpectrumCollection
from .utils import computation_wrapper


__all__ = ['moment']


def moment(spectrum, regions=None, order=0, axis=-1):
    """
    Estimate the moment of the spectrum.


    Parameters
    ----------
    spectrum : `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object over which the width will be calculated.

    regions: `~specutils.SpectralRegion` or list of `~specutils.SpectralRegion`
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
    if isinstance(spectrum, SpectrumCollection):
        return [computation_wrapper(_compute_moment, spec, regions,order=order, axis=axis)
                for spec in spectrum]
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

    dx = np.abs(np.diff(spectral_axis.bin_edges))
    m0 = np.sum(flux * dx, axis=axis)
    if order == 0:
        return m0

    dispersion = spectral_axis
    if len(flux.shape) > len(spectral_axis.shape):
        _shape = flux.shape[:-1] + (1,)
        dispersion = np.tile(spectral_axis, _shape)

    if order == 1:
        return np.sum(flux * dispersion * dx, axis=axis) / m0

    if order > 1:

        # By setting keepdims to True, the axes which are reduced are
        # left in the result as dimensions with size one. This means
        # that we can broadcast m1 correctly against dispersion.
        m1 = (np.sum(flux * dispersion * dx, axis=axis, keepdims=True)
              / np.sum(flux * dx, axis=axis, keepdims=True))

        return np.sum(flux * dx * (dispersion - m1) ** order, axis=axis) / m0
