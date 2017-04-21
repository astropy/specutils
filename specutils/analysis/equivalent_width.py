import numpy as np


__all__ = ['equivalent_width']

def equivalent_width(spectrum):
    """
    Does a naive equivalent width measures on the spectrum object.

    Parameters
    ----------
    spectrum : Spectrum1D
        The spectrum object overwhich the equivalent width will be calculated.

    Returns
    -------
    ew : float
        Equivalent width calculation.
    """
    # Continuum is always assumed to be 1.0
    avg_cont = np.median(spectrum.flux)

    # Average dispersion in the line region.
    avg_dx = np.mean(spectrum.dispersion[1:] - spectrum.dispersion[:-1])

    # Calculate equivalent width
    ew = ((avg_cont - spectrum.flux) * (avg_dx / avg_cont)).sum()

    return ew