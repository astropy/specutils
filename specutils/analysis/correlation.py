import astropy.units as u
import numpy as np
from astropy import constants as const
from astropy.nddata import StdDevUncertainty
from astropy.units import Quantity
from scipy.signal.windows import tukey
from scipy.signal import correlate

from ..manipulation import LinearInterpolatedResampler
from .. import Spectrum1D

__all__ = ['template_correlate', 'template_logwl_resample']

_KMS = u.Unit('km/s')  # for use below without having to create a composite unit


def template_correlate(observed_spectrum, template_spectrum, lag_units=_KMS,
                       apodization_window=0.5, resample=True, method="direct"):
    """
    Compute cross-correlation of the observed and template spectra.


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
    lag_units: `~astropy.units.Unit`
        Must be a unit with velocity physical type for lags in velocity. To
        output the lags in redshift, use ``u.dimensionless_unscaled``.
    apodization_window: float, callable, or None
        If a callable, will be treated as a window function for apodization of
        the cross-correlation (should behave like a `~scipy.signal.windows`
        window function, with ``sym=True``). If a float, will be treated as the
        ``alpha`` parameter for a Tukey window (`~scipy.signal.windows.tukey`),
        in units of pixels. If None, no apodization will be performed
    resample: bool or dict
        If True or a dictionary, resamples the spectrum and template following
        the process in `template_logwl_resample`. If a dictionary, it will be
        used as the keywords for `template_logwl_resample`.  For example,
        ``resample=dict(delta_log_wavelength=.1)`` would be the same as calling
        ``template_logwl_resample(spectrum, template, delta_log_wavelength=.1)``.
        If False, *no* resampling is performed (and the user is responsible for
        a sensible resampling).
    method: str
        If you choose "FFT", the correlation will be done through the use
        of convolution and will be calculated faster (for small spectral
        resolutions it is often correct), otherwise the correlation is determined
        directly from sums (the "direct" method in `~scipy.signal.correlate`).

    Returns
    -------
    (`~astropy.units.Quantity`, `~astropy.units.Quantity`)
        Arrays with correlation values and lags in km/s
    """

    # resample if the user requested to log wavelength
    if resample:
        if resample is True:
            resample_kwargs = dict()  # use defaults
        else:
            resample_kwargs = resample
        log_spectrum, log_template = template_logwl_resample(observed_spectrum,
                                                             template_spectrum,
                                                             **resample_kwargs)
    else:
        log_spectrum = observed_spectrum
        log_template = template_spectrum

    # apodize (might be a no-op if apodization_window is None)
    observed_log_spectrum, template_log_spectrum = _apodize(log_spectrum,
                                                            log_template,
                                                            apodization_window)
    # Normalize template
    normalization = _normalize(observed_log_spectrum, template_log_spectrum)

    # Not sure if we need to actually normalize the template. Depending
    # on the specific data uncertainty, the normalization factor
    # may turn out negative. That causes a flip of the correlation function,
    # in which the maximum (correlation peak) is no longer meaningful.
    if normalization < 0.:
        normalization = 1.

    # Correlate
    if method.lower() == "fft":
        corr = correlate(observed_log_spectrum.flux.value,
                         (template_log_spectrum.flux.value * normalization),
                         method="fft")
    else:
        corr = correlate(observed_log_spectrum.flux.value,
                         (template_log_spectrum.flux.value * normalization),
                         method="direct")

    # Compute lag
    # wave_l is the wavelength array equally spaced in log space.
    wave_l = observed_log_spectrum.spectral_axis.value
    delta_log_wave = np.log10(wave_l[1]) - np.log10(wave_l[0])
    deltas = (np.array(range(len(corr))) - len(corr)/2 + 0.5) * delta_log_wave
    lags = np.power(10., deltas) - 1.

    if u.dimensionless_unscaled.is_equivalent(lag_units):
        lags = Quantity(lags, u.dimensionless_unscaled)
    elif _KMS.is_equivalent(lag_units):
        lags = lags * const.c.to(lag_units)
    else:
        raise u.UnitsError('lag_units must be either velocity or dimensionless')

    return corr * u.dimensionless_unscaled, lags


def _apodize(spectrum, template, apodization_window):
    # Apodization. Must be performed after resampling.
    if apodization_window is None:
        clean_spectrum = spectrum
        clean_template = template
    else:
        if callable(apodization_window):
            window = apodization_window
        else:
            def window(wlen):
                return tukey(wlen, alpha=apodization_window)
        clean_spectrum = spectrum * window(len(spectrum.spectral_axis))
        clean_template = template * window(len(template.spectral_axis))

    return clean_spectrum, clean_template


def template_logwl_resample(spectrum, template, wblue=None, wred=None,
                            delta_log_wavelength=None,
                            resampler=LinearInterpolatedResampler()):
    """
    Resample a spectrum and template onto a common log-spaced spectral grid.

    If wavelength limits are not provided, the function will use
    the limits of the merged (observed+template) wavelength scale
    for building the log-wavelength scale.

    For the wavelength step, the function uses either the smallest wavelength
    interval found in the *observed* spectrum, or takes it from the
    ``delta_log_wavelength`` parameter.

    Parameters
    ----------
    observed_spectrum : :class:`~specutils.Spectrum1D`
        The observed spectrum.
    template_spectrum : :class:`~specutils.Spectrum1D`
        The template spectrum.
    wblue, wred: float
        Wavelength limits to include in the correlation.
    delta_log_wavelength: float
        Log-wavelength step to use to build the log-wavelength
        scale. If None, use limits defined as explained above.
    resampler
        A specutils resampler to use to actually do the resampling.  Defaults to
        using a `~specutils.manipulation.LinearInterpolatedResampler`.

    Returns
    -------
    resampled_observed : :class:`~specutils.Spectrum1D`
        The observed spectrum resampled to a common spectral_axis.
    resampled_template: :class:`~specutils.Spectrum1D`
        The template spectrum resampled to a common spectral_axis.
    """

    # Build an equally-spaced log-wavelength array based on
    # the input and template spectrum's limit wavelengths and
    # smallest sampling interval. Consider only the observed spectrum's
    # sampling, since it's the one that counts for the final accuracy
    # of the correlation. Alternatively, use the wred and wblue limits,
    # and delta log wave provided by the user.
    #
    # We work with separate float and units entities instead of Quantity
    # instances, due to the profusion of log10 and power function calls
    # (they only work on floats)
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

    return clean_spectrum, clean_template


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
    if hasattr(observed_spectrum, "uncertainty") and observed_spectrum.uncertainty is not None:
        unc = observed_spectrum.uncertainty.represent_as(StdDevUncertainty)
        num = np.nansum((observed_spectrum.flux*template_spectrum.flux)/(unc.array**2))
        denom = np.nansum((template_spectrum.flux/unc.array)**2)
    else:
        num = np.nansum(observed_spectrum.flux*template_spectrum.flux)
        denom = np.nansum(template_spectrum.flux**2)

    return num/denom
