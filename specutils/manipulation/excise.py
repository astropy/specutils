import numpy as np
import astropy.units as u

__all__ = ['excise_region']

def excise_region(spectrum, region):
    """
    Excise a spectral region from the spectrum.

    Parameters
    ----------
    spectrum: ``specutils.spectra.spectrum1d.Spectrum1D``
        The spectrum object from which the region will be excised.

    region: ``specutils.utils.SpectralRegion``
        The lower and upper disperion values (with ``astropy.units`` units)

    Return
    ------
    spectrum: ``specutils.spectra.spectrum1d.Spectrum1D``
        Excised spectrum.

    Notes
    -----
    The region excised is a discrete subset of the input spectrum. No interpolation is done
    on the left and right side of the spectrum.

    """

    dispersion = spectrum.spectral_axis

    left_index = np.nonzero(dispersion >= region.lower)
    if len(left_index) == 0:
        raise ValueError('No dispersion values greater than {}'.format(region.lower))
    else:
        left_index = left_index[0][0]

    right_index = np.nonzero(dispersion < region.upper)
    if len(right_index) == 0:
        raise ValueError('No dispersion values less than {}'.format(region.upper))
    else:
        right_index = right_index[0][-1]

    if left_index >= right_index:
        raise ValueError('Lower region, {}, appears to be greater than the upper region, {}.'.format(region[0], region[1]))

    return spectrum[left_index:right_index]

