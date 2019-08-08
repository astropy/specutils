import astropy.units as u
import numpy as np
from astropy.nddata import StdDevUncertainty

from ..spectra.spectrum1d import Spectrum1D
from ..spectra.spectrum_collection import SpectrumCollection
from ..analysis import template_match
from ..manipulation import FluxConservingResampler


def _normalize(observed_spectrum, template_spectrum):
    num = np.sum((observed_spectrum.flux*template_spectrum.flux)/(observed_spectrum.uncertainty.array**2))
    denom = np.sum((template_spectrum.flux/observed_spectrum.uncertainty.array)**2)
    normalized = num/denom

    return normalized

def get_chi_square(observed_spectrum, template_spectrum):
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

    return chi2

def test_template_match_spectrum():
    """
    Test template_match when both observed and template spectra have the same wavelength axis
    """
    # Create test spectra
    spec_axis = np.linspace(0, 50, 50) * u.AA
    spec = Spectrum1D(spectral_axis=spec_axis,
                      flux=np.random.randn(50) * u.Jy,
                      uncertainty=StdDevUncertainty(np.random.sample(50), unit='Jy'))

    spec1 = Spectrum1D(spectral_axis=spec_axis,
                       flux=np.random.randn(50) * u.Jy,
                       uncertainty=StdDevUncertainty(np.random.sample(50), unit='Jy'))

    # Get result from template_match
    tm_result = template_match.template_match(spec, spec1)

    # Calculate result independently of template_match
    result = get_chi_square(spec, spec1)

    assert result == tm_result[1]

def test_template_match_with_resample():
    """
    Test template_match when both observed and template spectra have different wavelength axis using resampling
    """
    # Create test spectra
    spec_axis1 = np.linspace(0, 50, 50) * u.AA
    spec_axis2 = np.linspace(0, 50, 50) * u.AA
    spec = Spectrum1D(spectral_axis=spec_axis1,
                      flux=np.random.randn(50) * u.Jy,
                      uncertainty=StdDevUncertainty(np.random.sample(50), unit='Jy'))

    spec1 = Spectrum1D(spectral_axis=spec_axis2,
                       flux=np.random.randn(50) * u.Jy,
                       uncertainty=StdDevUncertainty(np.random.sample(50), unit='Jy'))

    # Get result from template_match
    tm_result = template_match.template_match(spec, spec1)

    result = get_chi_square(spec, spec1)

    assert result == tm_result[1]

def test_template_match_list():
    """
    Test template_match when template spectra are in a list
    """
    # Create test spectra
    spec_axis1 = np.linspace(0, 50, 50) * u.AA
    spec_axis2 = np.linspace(0, 50, 50) * u.AA
    spec = Spectrum1D(spectral_axis=spec_axis1,
                      flux=np.random.randn(50) * u.Jy,
                      uncertainty=StdDevUncertainty(np.random.sample(50), unit='Jy'))

    spec1 = Spectrum1D(spectral_axis=spec_axis2,
                       flux=np.random.randn(50) * u.Jy,
                       uncertainty=StdDevUncertainty(np.random.sample(50), unit='Jy'))
    spec2 = Spectrum1D(spectral_axis=spec_axis2,
                       flux=np.random.randn(50) * u.Jy,
                       uncertainty=StdDevUncertainty(np.random.sample(50), unit='Jy'))

    # Combine spectra into list
    template_list = [spec1, spec2]

    # Get result from template_match
    tm_result = template_match.template_match(spec, template_list)

    # Loop through list in order to find out which template spectra has the smallest chi square
    chi2_min = None
    smallest_chi_spec = None
    for spectrum in template_list:
        chi2 = get_chi_square(spec, spectrum)
        if chi2_min is None or chi2 < chi2_min:
            chi2_min = chi2
            smallest_chi_spec = spectrum


    assert chi2_min == tm_result[1]
    assert smallest_chi_spec == tm_result[0]

def test_template_match_spectrum_collection():
    """
    Test template_match when template spectra are in a SpectrumCollection object
    """
    # Create test spectra
    spec_axis1 = np.linspace(0, 50, 50) * u.AA
    spec_axis2 = np.linspace(0, 50, 50) * u.AA
    spec = Spectrum1D(spectral_axis=spec_axis1,
                      flux=np.random.randn(50) * u.Jy,
                      uncertainty=StdDevUncertainty(np.random.sample(50), unit='Jy'))

    spec1 = Spectrum1D(spectral_axis=spec_axis2,
                       flux=np.random.randn(50) * u.Jy,
                       uncertainty=StdDevUncertainty(np.random.sample(50), unit='Jy'))
    spec2 = Spectrum1D(spectral_axis=spec_axis2,
                       flux=np.random.randn(50) * u.Jy,
                       uncertainty=StdDevUncertainty(np.random.sample(50), unit='Jy'))

    # Combine spectra into SpectrumCollection object
    spec_coll = SpectrumCollection.from_spectra([spec1, spec2])

    # Get result from template_match
    tm_result = template_match.template_match(spec, spec_coll)

    # Loop through SpectrumCollection in order to find out which template spectra has the smallest chi square
    chi2_min = None
    for spectrum in spec_coll:
        chi2 = get_chi_square(spec, spectrum)
        if chi2_min is None or chi2 < chi2_min:
            chi2_min = chi2

    assert chi2_min == tm_result[1]
