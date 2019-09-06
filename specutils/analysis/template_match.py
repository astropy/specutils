from ..spectra.spectrum1d import Spectrum1D
from ..spectra.spectrum_collection import SpectrumCollection
from ..manipulation import FluxConservingResampler
import numpy as np

def _normalize(observed_spectrum, template_spectrum):
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

def _template_match(observed_spectrum, template_spectrum):
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
    fluxc_resample = FluxConservingResampler()
    template_obswavelength = fluxc_resample(template_spectrum, observed_spectrum.wavelength)

    # Normalize spectra
    normalization = _normalize(observed_spectrum, template_obswavelength)

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

def template_match(observed_spectrum, spectral_templates):
    """
    Find what instance spectral_templates is and run _template_match accordingly. If two template_spectra have the same
    chi2, the first template is returned

    Parameters
    ----------
    observed_spectrum : :class:`~specutils.Spectrum1D`
        the observed spectrum
    spectral_templates : :class:`~specutils.Spectrum1D` or :class:`~specutils.SpectrumCollection` or `list`
        the template spectra, which will be resampled and normalized and compared to the observed_spectrum, where the
        smallest chi2 and normalized_template_spectrum will be returned

    Returns
    -------
    normalized_template_spectrum : :class:`~specutils.Spectrum1D`
        the template_spectrum that has been normalized
    chi2 : `float`
        the chi2 of the flux of the observed_spectrum and the flux of the normalized_template_spectrum
    """
    if isinstance(spectral_templates, Spectrum1D):
        normalized_spectral_template, chi2 = _template_match(observed_spectrum, spectral_templates)
        return normalized_spectral_template, chi2

    # Loop through spectra in list and return spectrum with lowest chi square
    # and its corresponding chi square
    elif isinstance(spectral_templates, list) or isinstance(spectral_templates, SpectrumCollection):
        chi2_min = None
        smallest_chi_spec = None

        for spectrum in spectral_templates:
            normalized_spectral_template, chi2 = _template_match(observed_spectrum, spectrum)
            if chi2_min is None or chi2 < chi2_min:
                chi2_min = chi2
                smallest_chi_spec = normalized_spectral_template

        return smallest_chi_spec, chi2_min
