import numpy as np
from astropy.nddata import StdDevUncertainty

from ..manipulation import (FluxConservingResampler,
                            LinearInterpolatedResampler,
                            SplineInterpolatedResampler)
from ..spectra.spectrum1d import Spectrum1D

__all__ = ['template_match', 'template_redshift']


def _normalize_for_template_matching(observed_spectrum, template_spectrum, stddev=None):
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
    if stddev is None:
        stddev = observed_spectrum.uncertainty.represent_as(StdDevUncertainty).quantity
    num = np.sum((observed_spectrum.flux*template_spectrum.flux) / (stddev**2))
    denom = np.sum((template_spectrum.flux / stddev)**2)

    return num/denom


def _resample(resample_method):
    """
    Find the user preferred method of resampling the template spectrum to fit
    the observed spectrum.

    Parameters
    ----------
    resample_method: `string`
        The type of resampling to be done on the template spectrum.

    Returns
    -------
    :class:`~specutils.ResamplerBase`
        This is the actual class that will handle the resampling.
    """
    if resample_method == "flux_conserving":
        return FluxConservingResampler()

    if resample_method == "linear_interpolated":
        return LinearInterpolatedResampler()

    if resample_method == "spline_interpolated":
        return SplineInterpolatedResampler()

    return None


def _chi_square_for_templates(observed_spectrum, template_spectrum, resample_method):
    """
    Resample the template spectrum to match the wavelength of the observed
    spectrum. Then, calculate chi2 on the flux of the two spectra.

    Parameters
    ----------
    observed_spectrum : :class:`~specutils.Spectrum1D`
        The observed spectrum.
    template_spectrum : :class:`~specutils.Spectrum1D`
        The template spectrum, which will be resampled to match the wavelength
        of the observed spectrum.

    Returns
    -------
    normalized_template_spectrum : :class:`~specutils.Spectrum1D`
        The normalized spectrum template.
    chi2 : `float`
        The chi2 of the flux of the observed spectrum and the flux of the
        normalized template spectrum.
    """
    # Resample template
    if _resample(resample_method) != 0:
        fluxc_resample = _resample(resample_method)
        template_obswavelength = fluxc_resample(template_spectrum,
                                                observed_spectrum.spectral_axis)

    # Convert the uncertainty to standard deviation if needed
    stddev = observed_spectrum.uncertainty.represent_as(StdDevUncertainty).array

    # Normalize spectra
    normalization = _normalize_for_template_matching(observed_spectrum,
                                                     template_obswavelength,
                                                     stddev)

    # Numerator
    num_right = normalization * template_obswavelength.flux
    num = observed_spectrum.flux - num_right

    # Denominator
    denom = stddev * observed_spectrum.flux.unit

    # Get chi square
    result = (num/denom)**2
    chi2 = np.sum(result.value)

    # Create normalized template spectrum, which will be returned with
    # corresponding chi2
    normalized_template_spectrum = Spectrum1D(
        spectral_axis=template_spectrum.spectral_axis,
        flux=template_spectrum.flux*normalization)

    return normalized_template_spectrum, chi2


def template_match(observed_spectrum, spectral_templates,
                   resample_method="flux_conserving",
                   redshift=None):
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
    resample_method : `string`
        Three resample options: flux_conserving, linear_interpolated, and spline_interpolated.
        Anything else does not resample the spectrum.
    redshift : 'float', `int`, `list`, `tuple`, 'numpy.array`
        If the user knows the redshift they want to apply to the spectrum/spectra within spectral_templates,
        then this float or int value redshift can be applied to each template before attempting the match.
        Or, alternatively, an iterable with redshift values to be applied to each template, before computation
        of the corresponding chi2 value, can be passed via this same parameter. For each template, the redshift
        value that results in the smallest chi2 is used.

    Returns
    -------
    normalized_template_spectrum : :class:`~specutils.Spectrum1D`
        The template spectrum that has been normalized.
    redshift : `None` or `float`
        The value of the redshift that provides the smallest chi2.
    smallest_chi_index : `int`
        The index of the spectrum with the smallest chi2 in spectral templates.
    chi2_min : `float`
        The chi2 of the flux of the observed_spectrum and the flux of the
        normalized template spectrum.
    chi2_list : `list`
        A list with all chi2 values found for each template spectrum.
    """
    final_redshift = None

    if hasattr(spectral_templates, 'flux') and len(spectral_templates.flux.shape) == 1:

        # Account for redshift if provided
        chi2_list = []
        if redshift is not None:
            results = template_redshift(observed_spectrum, spectral_templates, redshift=redshift)
            redshifted_spectrum, final_redshift, normalized_spectral_template, chi2, chi2_inner_list = results # noqa

        else:
            normalized_spectral_template, chi2 = _chi_square_for_templates(
                observed_spectrum, spectral_templates, resample_method)

        chi2_list.append(chi2)

        return normalized_spectral_template, final_redshift, 0, chi2, chi2_list

    # At this point, the template spectrum is either a ``SpectrumCollection``
    # or a multi-dimensional``Spectrum1D``. Loop through the object and return
    # the template spectrum with the lowest chi square and its corresponding
    # chi square.
    chi2_min = None
    smallest_chi_spec = None
    chi2_list = []

    final_redshift = None
    temp_redshift = None

    for index, spectrum in enumerate(spectral_templates):

        # Account for redshift if provided
        if redshift is not None:
            results = template_redshift(observed_spectrum, spectrum, redshift=redshift)
            redshifted_spectrum, temp_redshift, normalized_spectral_template, chi2, chi2_inner_list = results # noqa

            chi2_list.append(chi2_inner_list)

        else:
            normalized_spectral_template, chi2 = _chi_square_for_templates(
                observed_spectrum, spectrum, resample_method)

        if chi2_min is None or chi2 < chi2_min:
            chi2_min = chi2
            smallest_chi_spec = normalized_spectral_template
            smallest_chi_index = index

            final_redshift = temp_redshift

    return smallest_chi_spec, final_redshift, smallest_chi_index, chi2_min, chi2_list


def template_redshift(observed_spectrum, template_spectrum, redshift):
    """
    Find the best-fit redshift for template_spectrum to match observed_spectrum using chi2.

    Parameters
    ----------
    observed_spectrum : :class:`~specutils.Spectrum1D`
        The observed spectrum.
    template_spectrum : :class:`~specutils.Spectrum1D`
        The template spectrum, which will have it's redshift calculated.
    redshift : `float`, `int`, `list`, `tuple`, 'numpy.array`
        A scalar or iterable with the redshift values to test.

    Returns
    -------
    redshifted_spectrum: :class:`~specutils.Spectrum1D`
        A new Spectrum1D object which incorporates the template_spectrum with a spectral_axis
        that has been redshifted using the final_redshift.
    final_redshift : `float`
        The best-fit redshift for template_spectrum to match the observed_spectrum.
    normalized_template_spectrum : :class:`~specutils.Spectrum1D`
        The normalized spectrum template.
    chi2_min: `float`
        The smallest chi2 value that was found.
    chi2_list : `list`
        A list with the chi2 values corresponding to each input redshift value.
    """
    chi2_min = None
    final_redshift = None
    chi2_list = []

    redshift = np.array(redshift).reshape((np.array(redshift).size,))

    # Loop which goes through available redshift values and finds the smallest chi2
    for rs in redshift:

        # Create new redshifted spectrum and run it through the chi2 method
        redshifted_spectrum = Spectrum1D(spectral_axis=template_spectrum.spectral_axis*(1+rs),
                                         flux=template_spectrum.flux,
                                         uncertainty=template_spectrum.uncertainty,
                                         meta=template_spectrum.meta)

        normalized_spectral_template, chi2 = _chi_square_for_templates(
            observed_spectrum, redshifted_spectrum, "flux_conserving")

        chi2_list.append(chi2)

        # Set new chi2_min if suitable replacement is found
        if not np.isnan(chi2) and (chi2_min is None or chi2 < chi2_min):
            chi2_min = chi2
            final_redshift = rs

    return redshifted_spectrum, final_redshift, normalized_spectral_template, chi2_min, chi2_list
