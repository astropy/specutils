import numpy as np

import astropy.units as u
from astropy import constants as const


from ..manipulation import FluxConservingResampler


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
    Compute cross-correlation of the observed and template spectra.
    This assumes that both are resampled into a common wavelength
    scale before being submitted to this function.

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

    # Not sure if we need to actually normalize the template. Depending
    # on the specific data uncertainty, the normalization factor
    # may turn out negative. That causes a flip of the correlation function,
    # in which the maximum (correlation peak) is no longer meaningful.
    if normalization < 0.:
        normalization = 1.

    # Resample both spectrum and template to log wavelength. This is assuming
    # that both spectrum and template share a common spectral axis.
    resampler = FluxConservingResampler()
    observed_log_spectrum = _log_resampler(observed_spectrum, resampler)
    template_log_spectrum = _log_resampler(template_spectrum, resampler)

    # Correlate
    corr = np.correlate(observed_log_spectrum.flux.value,
                        (template_log_spectrum.flux.value * normalization),
                        mode='full')

    # Compute lag in km/s
    # wave_l is the wavelength array equally spaced in log space.
    wave_l = observed_log_spectrum.spectral_axis.value
    delta_log_wave = np.log10(wave_l[1]) - np.log10(wave_l[0])
    deltas = (np.array(range(len(corr))) - len(corr)/2 + 0.5) * delta_log_wave
    lags = (np.power(10., deltas) - 1. ) * const.c.to('km/s')

    return (corr * u.dimensionless_unscaled, lags)


def _log_resampler(spectrum, resampler):

    # This is coded independently of the fact that spectrum and template
    # are required to be expressed on a common wavelength axis. Depending
    # on how this work evolves, we can eventually simplify this code to
    # take advantage of the assumption of common wavelengths.

    # Build an equally-spaced log-wavelength array based on
    # the input spectrum's limit wavelengths and number of
    # samples.
    nsamples = len(spectrum.spectral_axis)
    w0 = np.log10(spectrum.spectral_axis[0].value)
    w1 = np.log10(spectrum.spectral_axis[-1].value)
    dw = (w1 - w0) / nsamples
    log_wave_array = np.ones(nsamples) * w0
    for i in range(nsamples):
        log_wave_array[i] += dw * i

    # Build it's corresponding wavelength array
    wave_array = np.power(10., log_wave_array) * spectrum.spectral_axis.unit

    # Resample spectrum into wavelength array so built
    resampled_spectrum = resampler(spectrum, wave_array)

    # Velocity convention and rest value are not preserved by resampler.
    # This might be a bug to squash? Note that this is in here for
    # formal consistency only, since the correlation code in this
    # module doesn't require these to be defined.
    resampled_spectrum = resampled_spectrum.with_velocity_convention(
        spectrum.velocity_convention)
    resampled_spectrum.rest_value = spectrum.rest_value

    return resampled_spectrum
