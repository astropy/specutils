import astropy.units as u
import numpy as np
from astropy import constants as const
from astropy.units import Quantity
from scipy.signal.windows import tukey

from ..manipulation import LinearInterpolatedResampler
from ..spectra.spectrum1d import Spectrum1D

__all__ = ['template_correlate']


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


def template_correlate(observed_spectrum, template_spectrum,
                       wblue=None, wred=None, delta_log_wavelength=None,
                       alpha=0.5, resample=True, lag_units=u.Unit('km/s')):
    """
    Compute cross-correlation of the observed and template spectra.

    Both observed and template spectra are re-sampled into a
    log-wavelength scale, unless the 'resample' keyword is set
    to False. In which case, they are assumed to have been previously
    resampled into a common log-wavelength scale (they share the same
    or identical spectral_axis attributes).

    If wavelength limits are not provided, the function will use
    the limits of the merged (observed+template) wavelength scale
    for building the log-wavelength scale.

    As for the wavelength step, the function uses either the smallest
    wavelength interval found in the observed spectrum, or takes
    it from the `delta_log_wavelength` parameter.

    After re-sampling into log-wavelength, both observed and template
    spectra are apodized by a Tukey window in order to minimize edge
    and consequent non-periodicity effects and thus decrease
    high-frequency power in the correlation function. To turn off the
    apodization, use alpha=0.


    Parameters
    ----------
    observed_spectrum : :class:`~specutils.Spectrum1D`
        The observed spectrum.
    template_spectrum : :class:`~specutils.Spectrum1D`
        The template spectrum, which will be correlated with
        the observed spectrum.
    wblue, wred: float
        Wavelength limits to include in the correlation
    delta_log_wavelength: float
        Log-wavelength step to use to build the log-wavelength
        scale. If None, use limits defined as explained above.
    alhpa: float, default=20.
        Apodization parameter (`~scipy.signal.windows.tukey`)
    resample: Boolean
        Re-sample spectrum and template into log-wavelength?
    lag_units: `~astropy.units.Unit`
        Default is 'km/s'. To output the lags in redshift,
        use `u.dimensionless_unscaled` instead

    Returns
    -------
    (`~astropy.units.Quantity`, `~astropy.units.Quantity`)
        Arrays with correlation values and lags in km/s
    """

    # Preprocess data:
    # - resample to log wavelength
    # - apodize
    observed_log_spectrum, template_log_spectrum = _preprocess(observed_spectrum,
                                                               template_spectrum,
                                                               wblue, wred,
                                                               delta_log_wavelength,
                                                               resample, alpha)
    # Normalize template
    normalization = _normalize(observed_log_spectrum, template_log_spectrum)

    # Not sure if we need to actually normalize the template. Depending
    # on the specific data uncertainty, the normalization factor
    # may turn out negative. That causes a flip of the correlation function,
    # in which the maximum (correlation peak) is no longer meaningful.
    if normalization < 0.:
        normalization = 1.

    # Correlate
    corr = np.correlate(observed_log_spectrum.flux.value,
                        (template_log_spectrum.flux.value * normalization),
                        mode='full')

    # Compute lag
    # wave_l is the wavelength array equally spaced in log space.
    wave_l = observed_log_spectrum.spectral_axis.value
    delta_log_wave = np.log10(wave_l[1]) - np.log10(wave_l[0])
    deltas = (np.array(range(len(corr))) - len(corr)/2 + 0.5) * delta_log_wave
    lags = np.power(10., deltas) - 1.

    if lag_units == u.dimensionless_unscaled:
        lags = Quantity(lags, u.dimensionless_unscaled)
    else:
        lags = lags * const.c.to(lag_units)

    return (corr * u.dimensionless_unscaled, lags)


def _preprocess(spectrum, template, wblue, wred, delta_log_wavelength, resample, alpha):
    if resample:
        # Build an equally-spaced log-wavelength array based on
        # the input and template spectrum's limit wavelengths and
        # smallest sampling interval. Consider only the observed spectrum's
        # sampling, since it's the one that counts for the final accuracy
        # of the correlation. Alternatively, use the wred and wblue limits,
        # and delta log wave provided by the user.
        if wblue:
            w0 = np.log10(wblue)
        else:
            ws0 = np.log10(spectrum.spectral_axis[0].value)
            wt0 = np.log10(template.spectral_axis[0].value)
            w0 = min(ws0, wt0)

        if wred:
            w1 = np.log10(wred)
        else:
            ws1 = np.log10(spectrum.spectral_axis[-1].value)
            wt1 = np.log10(template.spectral_axis[-1].value)
            w1 = max(ws1, wt1)

        if delta_log_wavelength is None:
            ds = np.log10(spectrum.spectral_axis.value[1:]) - np.log10(spectrum.spectral_axis.value[:-1])
            dw = ds[np.argmin(ds)]
        else:
            dw = delta_log_wavelength

        nsamples = int((w1 - w0) / dw)

        log_wave_array = np.ones(nsamples) * w0
        for i in range(nsamples):
            log_wave_array[i] += dw * i

        # Build the corresponding wavelength array
        wave_array = np.power(10., log_wave_array) * spectrum.spectral_axis.unit

        # Resample spectrum and template into wavelength array so built
        resampler = LinearInterpolatedResampler()
        resampled_spectrum = resampler(spectrum, wave_array)
        resampled_template = resampler(template, wave_array)

        # Resampler leaves Nans on flux bins that aren't touched by it.
        # We replace with zeros. This has the net effect of zero-padding
        # the spectrum and/or template so they exactly match each other,
        # wavelengthwise.
        clean_spectrum_flux = np.nan_to_num(resampled_spectrum.flux.value) * resampled_spectrum.flux.unit
        clean_template_flux = np.nan_to_num(resampled_template.flux.value) * resampled_template.flux.unit

        clean_spectrum = Spectrum1D(spectral_axis=resampled_spectrum.spectral_axis,
                            flux=clean_spectrum_flux,
                            uncertainty=resampled_spectrum.uncertainty,
                            velocity_convention='optical',
                            rest_value=spectrum.rest_value)
        clean_template = Spectrum1D(spectral_axis=resampled_template.spectral_axis,
                            flux=clean_template_flux,
                            uncertainty=resampled_template.uncertainty,
                            velocity_convention='optical',
                            rest_value=template.rest_value)
    else:
        clean_spectrum = spectrum
        clean_template = template

    # Apodization. Must be performed after resampling.
    clean_spectrum = clean_spectrum * tukey(len(clean_spectrum.wavelength), alpha=alpha)
    clean_template = clean_template * tukey(len(clean_template.wavelength), alpha=alpha)

    return clean_spectrum, clean_template
