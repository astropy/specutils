import astropy.units as u
import numpy as np
from astropy.nddata import StdDevUncertainty

from ..spectra.spectrum1d import Spectrum1D
from ..spectra.spectrum_collection import SpectrumCollection
from ..analysis import template_match
from astropy.tests.helper import quantity_allclose


# TODO: Add some tests that are outliers: where the observed and template do not overlap (what happens?),
# TODO: where there is minimal overlap (1 point)?

def test_template_match_spectrum():
    """
    Test template_match when both observed and template spectra have the same wavelength axis
    """
    # Seed np.random so that results are consistent
    np.random.seed(42)

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

    # Create new spectrum for comparison
    spec_result = Spectrum1D(spectral_axis=spec_axis,
                       flux=spec1.flux * template_match._normalize(spec, spec1))

    assert quantity_allclose(tm_result[0].flux, spec_result.flux, atol=0.01*u.Jy)
    assert tm_result[1] == 40093.28353756253

def test_template_match_with_resample():
    """
    Test template_match when both observed and template spectra have different wavelength axis using resampling
    """
    np.random.seed(42)

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

    # Create new spectrum for comparison
    spec_result = Spectrum1D(spectral_axis=spec_axis1,
                       flux=spec1.flux * template_match._normalize(spec, spec1))

    assert quantity_allclose(tm_result[0].flux, spec_result.flux, atol=0.01*u.Jy)
    assert tm_result[1] == 40093.28353756253

def test_template_match_list():
    """
    Test template_match when template spectra are in a list
    """
    np.random.seed(42)

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

    assert tm_result[1] == 40093.28353756253

def test_template_match_spectrum_collection():
    """
    Test template_match when template spectra are in a SpectrumCollection object
    """
    np.random.seed(42)

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

    assert tm_result[1] == 40093.28353756253