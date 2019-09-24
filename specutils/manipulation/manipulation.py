"""
A module for analysis tools dealing with uncertainties or error analysis in
spectra.
"""

import copy
import numpy as np
import operator

__all__ = ['snr_threshold']


def snr_threshold(spectrum, value, op=operator.gt):
    """
    Calculate the mean S/N of the spectrum based on the flux and uncertainty
    in the spectrum. This will be calculated over the regions, if they
    are specified.

    Parameters
    ----------
    spectrum : `~specutils.Spectrum1D`, `~specutils.SpectrumCollection` or `~astropy.nddata.NDData`
        The spectrum object overwhich the S/N threshold will be calculated.

    value: ``float``
        Threshold value to be applied to flux / uncertainty.

    op: One of operator.gt, operator.ge, operator.lt, operator.le
        The mathenatical operator to apply for thresholding.

    Returns
    -------
    spectrum: `~specutils.Spectrum1D`
        Ouput object with ``spectrum.mask`` set based on threshold.

    Notes
    -----
    The input object will need to have the uncertainty defined in order for the SNR
    to be calculated.

    """

    if not hasattr(spectrum, 'uncertainty') or spectrum.uncertainty is None:
        raise Exception("S/N thresholding requires the uncertainty be defined.")

    if op not in [operator.gt, operator.ge, operator.lt, operator.le]:
        raise ValueError('Threshold operator must be greater-than (operator.gt), greater-than-or-equal (operator.ge),'+
                         'less-than (operator.lt) or less-than-or-equal (operator.le)')

    # Spectrum1D
    if hasattr(spectrum, 'flux'):
        data = spectrum.flux

    # NDData
    elif hasattr(spectrum, 'data'):
        data = spectrum.data * (spectrum.unit if spectrum.unit is not None else 1)
    else:
        raise ValueError('Could not find data attribute.')

    mask = op((data / (spectrum.uncertainty.array*spectrum.uncertainty.unit)), value)

    spectrum_out = copy.copy(spectrum)
    spectrum_out._mask = mask

    return spectrum_out
