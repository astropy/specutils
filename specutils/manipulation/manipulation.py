"""
A module for analysis tools dealing with uncertainties or error analysis in
spectra.
"""

import copy
import operator

from astropy.nddata import StdDevUncertainty

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

    op: One of operator.gt, operator.ge, operator.lt, operator.le or
        the str equivalent '>', '>=', '<', '<='
        The mathematical operator to apply for thresholding.

    Returns
    -------
    spectrum: `~specutils.Spectrum1D`
        Output object with ``spectrum.mask`` set based on threshold.

    Notes
    -----
    The input object will need to have the uncertainty defined in order for the SNR
    to be calculated.

    """

    # Setup the mapping
    operator_mapping = {
        '>': operator.gt,
        '<': operator.lt,
        '>=': operator.ge,
        '<=': operator.le
    }

    if not hasattr(spectrum, 'uncertainty') or spectrum.uncertainty is None:
        raise Exception("S/N thresholding requires the uncertainty be defined.")

    if (op not in [operator.gt, operator.ge, operator.lt, operator.le] and
            op not in operator_mapping.keys()):
        raise ValueError('Threshold operator must be a string or operator that represents ' +
                         'greater-than, less-than, greater-than-or-equal or ' +
                         'less-than-or-equal')

    # If the operator passed in is a string, then map to the
    # operator method.
    if isinstance(op, str):
        op = operator_mapping[op]

    # Spectrum1D
    if hasattr(spectrum, 'flux'):
        data = spectrum.flux

    # NDData
    elif hasattr(spectrum, 'data'):
        data = spectrum.data * (spectrum.unit if spectrum.unit is not None else 1)
    else:
        raise ValueError('Could not find data attribute.')

    # NDData convention: Masks should follow the numpy convention that valid
    # data points are marked by False and invalid ones with True.
    mask = ~op(data / (spectrum.uncertainty.represent_as(StdDevUncertainty).quantity), value)

    spectrum_out = copy.copy(spectrum)
    spectrum_out._mask = mask

    return spectrum_out
