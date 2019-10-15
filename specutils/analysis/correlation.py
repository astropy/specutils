import numpy as np

from ..manipulation import (FluxConservingResampler,
                            LinearInterpolatedResampler,
                            SplineInterpolatedResampler)
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
    # Normalize spectra
    normalization = _normalize(observed_spectrum, template_spectrum)

    # Numerator
    num_right = normalization * template_spectrum.flux
    num = observed_spectrum.flux - num_right

    # Denominator
    denom = observed_spectrum.uncertainty.array * observed_spectrum.flux.unit

    # Get chi square
    result = (num/denom)**2
    chi2 = np.sum(result.value)

    # Create normalized template spectrum, which will be returned with
    # corresponding chi2
    normalized_template_spectrum = Spectrum1D(
        spectral_axis=template_spectrum.spectral_axis,
        flux=template_spectrum.flux*normalization)

    return normalized_template_spectrum, chi2


# TODO not sure if this is part of the core functionality, or is
# something the user will do outside.

def correlation(observed_spectrum, spectral_templates):
    """
    Find which spectral templates is the best fit to an observed spectrum by
    computing the chi-squared. If two template_spectra have the same chi2, the
    first template is returned.

    Parameters
    ----------
    observed_spectrum : :class:`~specutils.Spectrum1D`
        The observed spectrum.
    spectral_templates : :class:`~specutils.Spectrum1D` or :class:`~specutils.SpectrumCollection` or `list`
        That will give a single :class:`~specutils.Spectrum1D` when iterated
        over. The template spectra, which will be resampled, normalized, and
        compared to the observed spectrum, where the smallest chi2 and
        normalized template spectrum will be returned.

    Returns
    -------
    normalized_template_spectrum : :class:`~specutils.Spectrum1D`
        The template spectrum that has been normalized.
    chi2 : `float`
        The chi2 of the flux of the observed_spectrum and the flux of the
        normalized template spectrum.
    smallest_chi_index : `int`
        The index of the spectrum with the smallest chi2 in spectral templates.
    """
    if hasattr(spectral_templates, 'flux') and len(spectral_templates.flux.shape) == 1:
        normalized_spectral_template, chi2 = template_correlate(
            observed_spectrum, spectral_templates)

        return normalized_spectral_template, chi2

    # At this point, the template spectrum is either a ``SpectrumCollection``
    # or a multi-dimensional``Spectrum1D``. Loop through the object and return
    # the template spectrum with the lowest chi square and its corresponding
    # chi square.
    chi2_min = None
    smallest_chi_spec = None

    for index, spectrum in enumerate(spectral_templates):
        normalized_spectral_template, chi2 = template_correlate(
            observed_spectrum, spectrum)

        if chi2_min is None or chi2 < chi2_min:
            chi2_min = chi2
            smallest_chi_spec = normalized_spectral_template
            smallest_chi_index = index

    return smallest_chi_spec, chi2_min, smallest_chi_index
