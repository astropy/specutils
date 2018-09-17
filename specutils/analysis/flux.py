from __future__ import division

import numpy as np
from .utils import computation_wrapper


__all__ = ['line_flux', 'equivalent_width']


def line_flux(spectrum, region=None):
    """
    Computes the line flux over a spectrum. Applies a region if given.

    Parameters
    ----------
    spectrum : Spectrum1D
        The spectrum object over which the line flux will be calculated.

    Returns
    -------
    flux : `~astropy.units.Quantity`
        line flux result
    """
    return computation_wrapper(_compute_line_flux, spectrum, region)


def equivalent_width(spectrum, region=None):
    """
    Does a naive equivalent width measures on the spectrum object.

    Parameters
    ----------
    spectrum : Spectrum1D
        The spectrum object overwhich the equivalent width will be calculated.

    Returns
    -------
    ew : `~astropy.units.Quantity`
        Equivalent width calculation.

    TODO:  what frame of reference do you want the spectral_axis to be in ???
    """
    return computation_wrapper(_compute_equivalent_width, spectrum, region)


def _compute_line_flux(spectrum, region=None):

    if region is not None:
        calc_spectrum = region.extract(spectrum)
    else:
        calc_spectrum = spectrum

    # Average dispersion in the line region
    avg_dx = np.diff(spectrum.spectral_axis)

    line_flux = np.sum(calc_spectrum.flux[1:] * avg_dx)

    # TODO: we may want to consider converting to erg / cm^2 / sec by default
    return line_flux


def _compute_equivalent_width(spectrum, region=None):

    if region is not None:
        calc_spectrum = region.extract(spectrum)
    else:
        calc_spectrum = spectrum

    spectral_axis = calc_spectrum.spectral_axis
    dx = spectral_axis[-1] - spectral_axis[0]

    # Calculate equivalent width
    ew =  dx - (_compute_line_flux(spectrum, region) / calc_spectrum.flux.unit)

    return ew
