import numpy as np

import astropy.units as u

from ..spectra.spectrum1d import Spectrum1D

def _normalize(observed_spectrum, template_spectrum):
    """
    Calculate a scale factor to be applied to the template spectrum so the
    total flux in both spectra will be the same.

    Parameters
    ----------
    observed_spectrum : :class:`~specutils.Spectrum1D`
        The observed spectrum.
    template_spectrum : :class:`~specutils.Spectrum1D`
        The template spectrum, which needs to be normalized in order to be
        compared with the observed spectrum.

    Returns
    -------
    `float`
        A float which will normalize the template spectrum's flux so that it
        can be compared to the observed spectrum.
    """
    num = np.sum((observed_spectrum.flux*template_spectrum.flux)/
                 (observed_spectrum.uncertainty.array**2))
    denom = np.sum((template_spectrum.flux/
                    observed_spectrum.uncertainty.array)**2)

    return num/denom


def template_correlate(observed_spectrum, template_spectrum):
    """
    Compute cross-correlation of the observed and template spectra

    Parameters
    ----------
    observed_spectrum : :class:`~specutils.Spectrum1D`
        The observed spectrum.
    template_spectrum : :class:`~specutils.Spectrum1D`
        The template spectrum, which will be correlated with
        the observed spectrum.

    Returns
    -------
    correlation function : :class:`~specutils.Spectrum1D`
        The normalized correlation function.
    """
    # Normalize template
    normalization = _normalize(observed_spectrum, template_spectrum)

    # Correlation
    corr = np.correlate(observed_spectrum.flux.value,
                        (template_spectrum.flux.value * normalization),
                        mode='full')

    # Retun correlation function as a Spectrum1D instance
    lags = np.array([(a-len(corr)/2+0.5) for a in range(len(corr))]) * \
           observed_spectrum.spectral_axis.unit

    spectrum_midpoint = int(len(observed_spectrum.spectral_axis) / 2)
    lags *= (observed_spectrum.spectral_axis[spectrum_midpoint].value -
             observed_spectrum.spectral_axis[spectrum_midpoint-1].value)

    correlation_function = Spectrum1D(spectral_axis=lags,
                                      flux=corr * u.dimensionless_unscaled)

    return correlation_function

