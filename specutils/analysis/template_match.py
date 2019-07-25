from ..spectra.spectrum1d import Spectrum1D
from ..spectra.spectrum_collection import SpectrumCollection
from ..manipulation import FluxConservingResample
from scipy.stats import chisquare
import numpy as np
import math

def _normalize(observed_spectrum, template_spectrum):
    print(observed_spectrum.uncertainty)

    normalization = math.pow(np.sum((observed_spectrum*template_spectrum)/(observed_spectrum.uncertainty)), 2) \
                    / math.pow(np.sum((template_spectrum/observed_spectrum.uncertainty)), 2)
    return normalization

def _template_match(observed_spectrum, template_spectrum):
    """
    Resample the template spectrum to match the wavelength of the observed spectrum.
    Then, run scipy.stats.chisquare on the flux of the two spectra.

    :param observed_spectrum: The observed spectrum
    :param template_spectrum: The template spectrum, which will be resampled to match the wavelength of the observed
        spectrum
    :return: chi square of the flux of the observed spectrum and the flux of the template spectrum
    """
    fluxc_resample = FluxConservingResample()
    template_spectrum1D = fluxc_resample(template_spectrum, observed_spectrum.wavelength)

    # calculate chi square
    # x2 = chisquare(observed_spectrum.flux, template_spectrum1D.flux)

    normalization = _normalize(observed_spectrum, template_spectrum1D)

    chi2 = np.sum(((observed_spectrum.flux-normalization*template_spectrum1D.flux)/observed_spectrum.uncertainty)**2 )
    return chi2

def template_match(observed_spectrum, spectral_templates):
    """
    Find what instance collection is and run _template_match accordingly
    """
    if isinstance(spectral_templates, Spectrum1D):
        return Spectrum1D, _template_match(observed_spectrum, spectral_templates)

    elif isinstance(spectral_templates, SpectrumCollection):
        pass

    # Loop through spectra in list and return spectrum with lowest chi square
    # and its corresponding chi square
    elif isinstance(spectral_templates, list):
        chi2_min = None
        smallest_chi_spec = None
        for spectrum in spectral_templates:
            chi2 = _template_match(observed_spectrum, spectrum)
            if chi2_min is None or chi2 < chi2_min:
                chi2_min = chi2
                smallest_chi_spec = spectrum

        return smallest_chi_spec, chi2_min
