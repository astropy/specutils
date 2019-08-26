from ..spectra.spectrum1d import Spectrum1D
from ..spectra.spectrum_collection import SpectrumCollection
from ..manipulation import FluxConservingResampler
import numpy as np

def _normalize(observed_spectrum, template_spectrum):
    """
    The normalization is necessary to bring the template to the same magnitude as the observation and minimize the chi^2
    """
    num = np.sum((observed_spectrum.flux*template_spectrum.flux)/(observed_spectrum.uncertainty.array**2))
    denom = np.sum((template_spectrum.flux/observed_spectrum.uncertainty.array)**2)
    normalized = num/denom

    return normalized

def _template_match(observed_spectrum, template_spectrum):
    """
    Resample the template spectrum to match the wavelength of the observed spectrum.
    Then, run chisquare on the flux of the two spectra.

    :param observed_spectrum: The observed spectrum
    :param template_spectrum: The template spectrum, which will be resampled to match the wavelength of the observed
        spectrum
    :return: chi square of the flux of the observed spectrum and the flux of the template spectrum
    """
    # Resample template
    fluxc_resample = FluxConservingResampler()
    template_spectrum1D = fluxc_resample(template_spectrum, observed_spectrum.wavelength)

    # Normalize spectra
    normalization = _normalize(observed_spectrum, template_spectrum1D)

    # Numerator
    num_right = normalization*template_spectrum1D.flux
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
    Find what instance collection is and run _template_match accordingly
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
