"""Functions for spectral statistical analysis."""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# THIRD-PARTY
import numpy as np

# LOCAL
from ..core.data import Data

__all__ = ['extract', 'stats', 'eq_width']


def extract(data, x_range):
    """Extract a region from a spectrum.

    Parameters
    ----------
    data : ``Data``
        Contains the spectrum to be extracted.

    x_range : tuple
        A spectral coordinate range as in ``(wave1, wave2)``.

    Returns
    -------
    result : ``Data``
        Spectrum data with extracted region.

    """
    x = data.data
    y = data.dispersion

    slice = (x >= x_range[0]) & (x < x_range[1])

    result = Data(x[slice])
    result.set_x(x[slice], unit=spectrum_data.x.unit)
    result.set_y(y[slice], unit=spectrum_data.y.unit)

    return result


def stats(data):
    """Compute basic statistics for a spectral region
    contained in a ``Data`` instance.

    Parameters
    ----------
    data : ``Data``
        Typically this is returned by the :func:`extract` function.

    Returns
    -------
    statistics : dict
        Statistics results.

    """
    return {'mean':    np.mean(data),
            'median':  np.median(data),
            'stddev':  np.std(data),
            'total':   np.trapz(data),
            'npoints': len(data)}


def eq_width(cont1_stats, cont2_stats, line):
    """Compute an equivalent width given stats for two continuum
    regions, and a ``Data`` instance with the extracted
    spectral line region.

    This uses for now a very simple continuum subtraction method; i.e.,
    it just subtracts a constant from the line spectrum, where the
    constant is ``(continuum1[mean] + continuum2[mean]) / 2``.

    Parameters
    ----------
    cont1_stats, cont2_stats : dict
        This is returned by the :func:`stats` function.

    line : ``Data``
        This is returned by the :func:`extract` function.

    Returns
    -------
    flux, ew : float
        Flux and equivalent width values.

    """
    # average of 2 continuum regions.
    avg_cont = (cont1_stats['mean'] + cont2_stats['mean']) / 2.0

    # average dispersion in the line region.
    avg_dx = np.mean(line.x.data[1:] - line.x.data[:-1])

    # flux
    flux = np.sum(line.y.data - avg_cont) * avg_dx

    #  EW = Sum( (Fc-Fl)/Fc * dw
    ew =  np.sum((avg_cont - line.y.data) / avg_cont * avg_dx)

    return flux, ew
