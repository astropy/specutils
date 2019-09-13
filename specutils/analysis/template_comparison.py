from ..spectra.spectrum1d import Spectrum1D
from ..spectra.spectrum_collection import SpectrumCollection
from ..manipulation import FluxConservingResampler
from ..manipulation import LinearInterpolatedResampler
from ..manipulation import SplineInterpolatedResampler
import numpy as np

def _normalize_for_template_matching(observed_spectrum, template_spectrum):
    """
    Calculate a scale factor to be applied to the template spectrum so the total flux in both spectra will be the same.

    Parameters
    ----------
    observed_spectrum : :class:`~specutils.Spectrum1D`
        the observed spectrum
    template_spectrum : :class:`~specutils.Spectrum1D`
        the template spectrum, which needs to be normalized in order to be compared with the observed_spectrum

    Returns
    -------
    `float`
        A float which will normalize the template spectrum's flux so that it can be compared to the observed spectrum
    """
    num = np.sum((observed_spectrum.flux*template_spectrum.flux)/(observed_spectrum.uncertainty.array**2))
    denom = np.sum((template_spectrum.flux/observed_spectrum.uncertainty.array)**2)
    return num/denom

def _resample(resample_method):
    """
    Find the user preferred method of resampling the template spectrum to fit the observed spectrum.

    Parameters
    ----------
    resample_method: `string`
        The type of resampling to be done on the template spectrum

    Returns
    -------
    :class:`~specutils.ResamplerBase`
        This is the actual class that will handle the resampling
    """
    if resample_method == "flux_conserving":
        return FluxConservingResampler()
    elif resample_method == "linear_interpolated":
        return LinearInterpolatedResampler()
    elif resample_method == "spline_interpolated":
        return SplineInterpolatedResampler()
    else:
        return None

def _template_match(observed_spectrum, template_spectrum, resample_method):
    """
    Resample the template spectrum to match the wavelength of the observed spectrum.
    Then, calculate chi2 on the flux of the two spectra.

    Parameters
    ----------
    observed_spectrum : :class:`~specutils.Spectrum1D`
        the observed spectrum
    template_spectrum : :class:`~specutils.Spectrum1D`
        the template spectrum, which will be resampled to match the wavelength of the observed_spectrum

    Returns
    -------
    normalized_template_spectrum : :class:`~specutils.Spectrum1D`
        the template_spectrum that has been normalized
    chi2 : `float`
        the chi2 of the flux of the observed_spectrum and the flux of the normalized_template_spectrum
    """
    # Resample template
    if _resample(resample_method) != 0:
        fluxc_resample = _resample(resample_method)
        template_obswavelength = fluxc_resample(template_spectrum, observed_spectrum.wavelength)

    # Normalize spectra
    normalization = _normalize_for_template_matching(observed_spectrum, template_obswavelength)

    # Numerator
    num_right = normalization*template_obswavelength.flux
    num = observed_spectrum.flux - num_right

    # Denominator
    denom = observed_spectrum.uncertainty.array * observed_spectrum.flux.unit

    # Get chi square
    result = (num/denom)**2
    chi2 = np.sum(result)

    # Create normalized template spectrum, which will be returned with corresponding chi2
    normalized_template_spectrum = Spectrum1D(spectral_axis=template_spectrum.spectral_axis,
                                              flux=template_spectrum.flux*normalization)

    return normalized_template_spectrum, chi2

def template_match(observed_spectrum, spectral_templates, resample_method="flux_conserving"):
    """
    Find which spectral templates is the best fit to an observed spectrum by computing the chi-squared. If two
    template_spectra have the same chi2, the first template is returned

    Parameters
    ----------
    observed_spectrum : :class:`~specutils.Spectrum1D`
        the observed spectrum
    spectral_templates : :class:`~specutils.Spectrum1D` or :class:`~specutils.SpectrumCollection` or `list` or anything
        that will give a single :class:`~specutils.Spectrum1D` when iterated over.
        The template spectra, which will be resampled and normalized and compared to the observed_spectrum, where the
        smallest chi2 and normalized_template_spectrum will be returned

    Returns
    -------
    normalized_template_spectrum : :class:`~specutils.Spectrum1D`
        the template_spectrum that has been normalized
    chi2 : `float`
        the chi2 of the flux of the observed_spectrum and the flux of the normalized_template_spectrum
    smallest_chi_index : `int`
        the index of the spectrum with the smallest chi2 in spectral_templates
    """
    if hasattr(spectral_templates, 'flux') and len(spectral_templates.flux.shape)==1:
        normalized_spectral_template, chi2 = _template_match(observed_spectrum, spectral_templates, resample_method)
        return normalized_spectral_template, chi2

    # Loop through spectra in list and return spectrum with lowest chi square
    # and its corresponding chi square
    else:
        chi2_min = None
        smallest_chi_spec = None

        index = 0
        try:
            for spectrum in spectral_templates:
                normalized_spectral_template, chi2 = _template_match(observed_spectrum, spectrum, resample_method)
                if chi2_min is None or chi2 < chi2_min:
                    chi2_min = chi2
                    smallest_chi_spec = normalized_spectral_template
                    smallest_chi_index = index
                index += 1
        except Exception as e:
            print("Parameter spectral_templates is not iterable. The following error was fired: {}".format(e))

        return smallest_chi_spec, chi2_min, smallest_chi_index
