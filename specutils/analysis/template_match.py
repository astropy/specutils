from ..spectra.spectrum1d import Spectrum1D
from ..spectra.spectrum_collection import SpectrumCollection
from ..manipulation import FluxConservingResample
from scipy.stats import chisquare

def _template_match(observed_spectrum, template_spectrum):
    # resample the template spectrum to match the wavelength of the observed spectrum
    fluxc_resample = FluxConservingResample()
    output_spectrum1D = fluxc_resample(template_spectrum, observed_spectrum.wavelength)

    # calculate chi square
    x2 = chisquare(observed_spectrum.flux, output_spectrum1D.flux)
    return x2

def template_match(observed_spectrum, collection):
    """
    Find what instance collection is and run _template_match accordingly
    """
    if isinstance(collection, Spectrum1D):
        return Spectrum1D, _template_match(observed_spectrum, collection)

    elif isinstance(collection, SpectrumCollection):
        pass

    # Loop through spectra in list and return spectrum with lowest chi square
    # and its corresponding chi square
    elif isinstance(collection, list):
        min = None
        smallest_chi_spec = None
        for spectrum in collection:
            chi = _template_match(observed_spectrum, spectrum)
            if min is None or chi < min:
                min = chi
                smallest_chi_spec = spectrum

        return smallest_chi_spec, min
