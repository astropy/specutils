"""Functions for spectral statistical analysis."""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# THIRD-PARTY
import numpy as np
import scipy as sp

# LOCAL
from ..core.data import Data

__all__ = ['extract', 'stats', 'eq_width', 'fwzi', 'centroid']


def extract(data, x_range):
    """Extract a region from a spectrum.

    Parameters
    ----------
    data : `~specviz.core.data.Data`
        Contains the spectrum to be extracted.

    x_range : tuple
        A spectral coordinate range as in ``(wave1, wave2)``.

    Returns
    -------
    result : `~specviz.core.data.Data`
        Spectrum data with extracted region.

    Examples
    --------
    >>> d = Data(...)
    >>> d2 = extract(d, (10000, 20000))

    """
    y = data.data
    x = data.dispersion

    slice = (x >= x_range[0]) & (x < x_range[1])

    result = Data(y[slice], dispersion=x[slice])

    # TODO: Need to implement unit handling in Data
    #result.set_x(x[slice], unit=spectrum_data.x.unit)
    #result.set_y(y[slice], unit=spectrum_data.y.unit)

    return result


def stats(data):
    """Compute basic statistics for a spectral region
    contained in a `~specviz.core.data.Data` instance.

    Parameters
    ----------
    data : `~specviz.core.data.Data`
        Typically this is returned by the :func:`extract` function.

    Returns
    -------
    statistics : dict
        Statistics results.

    Examples
    --------
    >>> d = Data(...)
    >>> d_stats = stats(d)

    """
    if isinstance(data, Data):
        y = data.data
    else:
        y = data

    return {'mean':    np.mean(y),
            'median':  np.median(y),
            'stddev':  np.std(y),
            'total':   np.trapz(y),
            'npoints': len(y)}


def eq_width(cont1_stats, cont2_stats, line):
    """Compute an equivalent width given stats for two continuum
    regions, and a `~specviz.core.data.Data` instance with the extracted
    spectral line region.

    This uses for now a very simple continuum subtraction method; i.e.,
    it just subtracts a constant from the line spectrum, where the
    constant is ``(continuum1[mean] + continuum2[mean]) / 2``.

    Parameters
    ----------
    cont1_stats, cont2_stats : dict
        This is returned by the :func:`stats` function.

    line : `~specviz.core.data.Data`
        This is returned by the :func:`extract` function.

    Returns
    -------
    flux, ew : float
        Flux and equivalent width values.

    Examples
    --------
    >>> d = Data(...)
    >>> cont1 = extract(d, (100, 5000))
    >>> cont2 = extract(d, (18000, 20000))
    >>> cont1_stats = stats(cont1)
    >>> cont2_stats = stats(cont2)
    >>> line = extract(d, (15000, 18000))
    >>> flux, ew = eq_width(cont1_stats, cont2_stats, line)

    """
    # average of 2 continuum regions.
    avg_cont = (cont1_stats['mean'] + cont2_stats['mean']) / 2.0

    # average dispersion in the line region.
    avg_dx = np.mean(line.dispersion[1:] - line.dispersion[:-1])

    # flux
    flux = np.sum(line.data - avg_cont) * avg_dx

    #  EW = Sum( (Fc-Fl)/Fc * dw
    ew = np.abs(np.sum((avg_cont - line.data) / avg_cont * avg_dx))

    return flux, ew


# TODO: Can this be improved?
def fwzi(cont1_stats, cont2_stats, line):
    """Compute full width at zero intensity (FWZI) for the given spectrum.
    Continuum calculations are similar to :func:`eq_width`.

    Parameters
    ----------
    cont1_stats, cont2_stats : dict
        This is returned by the :func:`stats` function.

    line : `~specviz.core.data.Data`
        This is returned by the :func:`extract` function.

    Returns
    -------
    fwzi_value : float
        FWZI value.

    w_range : tuple
        Wavelengths used to calculate FWZI.

    Examples
    --------
    >>> d = Data(...)
    >>> cont1 = extract(d, (100, 5000))
    >>> cont2 = extract(d, (18000, 20000))
    >>> cont1_stats = stats(cont1)
    >>> cont2_stats = stats(cont2)
    >>> line = extract(d, (15000, 18000))
    >>> fwzi_value, w_range = fwzi(cont1_stats, cont2_stats, line)

    """
    # average of 2 continuum regions.
    avg_cont = (cont1_stats['mean'] + cont2_stats['mean']) * 0.5

    # Subtract continuum
    flux = line.data - avg_cont
    idx = len(flux) // 2

    # Guess absorption or emission -- dangerous?
    is_em = flux[idx] > avg_cont

    # Find wavelengths of zero intensity on each side
    w1 = line.dispersion[:idx]
    f1 = flux[:idx]
    w2 = line.dispersion[idx:]
    f2 = flux[idx:]
    if is_em:
        mask1 = f1 <= 0
        mask2 = f2 <= 0
    else:
        mask1 = f1 >= 0
        mask2 = f2 >= 0

    # If given line does not reach zero intensity, just take its edges
    if not np.any(mask1):
        mask1 = f1 == f1[0]
    if not np.any(mask2):
        mask2 = f2 == f2[-1]

    vmin = np.max(w1[mask1])
    vmax = np.min(w2[mask2])

    return vmax - vmin, (vmin, vmax)


def centroid(data):
    """Compute centroid for the given spectrum. ::

        w_cen = integral(wave*flux) / integral(flux)

    Parameters
    ----------
    data : `~specviz.core.data.Data`
        Extracted spectrum data.

    Returns
    -------
    wcen : float
        Centroid wavelength.

    Examples
    --------
    >>> d = Data(...)
    >>> line = extract(d, (15000, 20000))
    >>> wcen_em = centroid(line)

    """
    flux = data.data
    wave = data.dispersion

    wcen = sp.integrate.trapz(wave * flux, wave) / sp.integrate.trapz(flux,
                                                                      wave)

    return wcen
