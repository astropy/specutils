import numpy as np

import astropy.units as u


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
    num = np.nansum((observed_spectrum.flux*template_spectrum.flux)/
                 (observed_spectrum.uncertainty.array**2))
    denom = np.nansum((template_spectrum.flux/
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
    (`~astropy.units.Quantity`, `~astropy.units.Quantity`)
        Arrays with correlation values and lags in km/s
    """
    # Normalize template
    normalization = _normalize(observed_spectrum, template_spectrum)

    # Correlation
    corr = np.correlate(observed_spectrum.flux.value,
                        (template_spectrum.flux.value * normalization),
                        mode='full')

    # Lag in km/s
    equiv = getattr(u.equivalencies, 'doppler_{0}'.format(
        observed_spectrum.velocity_convention))(observed_spectrum.rest_value)

    lags = observed_spectrum.spectral_axis.to(u.km / u.s, equivalencies=equiv)
    diff = lags[1:] - lags[:-1]
    lags = (np.array(range(len(corr))) - len(corr)/2 + 0.5) * diff[0]

    return (corr * u.dimensionless_unscaled, lags)
