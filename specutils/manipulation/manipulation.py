"""
A module for analysis tools dealing with uncertainties or error analysis in
spectra.
"""

import copy
import numpy as np

__all__ = ['snr_threshold']


def snr_threshold(spectrum, value):
    """
    Calculate the mean S/N of the spectrum based on the flux and uncertainty
    in the spectrum. This will be calculated over the regions, if they
    are specified.

    Parameters
    ----------
    spectrum : `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object overwhich the equivalent width will be calculated.

    value: ``float``
        Threshold value to be applied to flux / uncertainty.

    Returns
    -------
    spectrum: `~specutils.Spectrum1D`
        Ouput spectrum with ``spectrum.mask`` set based on threshold.

    Notes
    -----
    The spectrum will need to have the uncertainty defined in order for the SNR
    to be calculated.

    """

    if not hasattr(spectrum, 'uncertainty') or spectrum.uncertainty is None:
        raise Exception("S/N thresholding requires the uncertainty be defined.")

    if hasattr(spectrum, 'flux'):
        data = spectrum.flux
    elif hasattr(spectrum, 'data'):
        data = spectrum.data * (spectrum.unit if spectrum.unit is not None else 1)
    else:
        raise ValueError('Could not find data attribute.')

    mask = (data / (spectrum.uncertainty.array*spectrum.uncertainty.unit)) > value

    spectrum_out = copy.copy(spectrum)
    spectrum_out._mask = mask

    return spectrum_out
