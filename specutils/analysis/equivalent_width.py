from __future__ import division

import numpy as np
from .utils import computation_wrapper


__all__ = ['equivalent_width']


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


def _compute_equivalent_width(spectrum, region=None):

    if region is not None:
        calc_spectrum = region.extract(spectrum)
    else:
        calc_spectrum = spectrum

    # Continuum is always assumed to be 1.0
    avg_cont = np.median(spectrum.flux)

    # Average dispersion in the line region
    avg_dx = np.mean(spectrum.wavelength[1:] - spectrum.wavelength[:-1])

    # Calculate equivalent width
    ew = ((avg_cont - spectrum.flux) * (avg_dx / avg_cont)).sum()

    return ew
