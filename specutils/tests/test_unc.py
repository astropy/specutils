import pytest
import numpy as np

import astropy.units as u
from astropy.tests.helper import quantity_allclose
from astropy.modeling import models
from astropy.nddata.nduncertainty import StdDevUncertainty, VarianceUncertainty, InverseVariance

from ..spectra import Spectrum1D, SpectralRegion
from ..manipulation import noise_region_uncertainty


def pure_noise_spectrum(amplitude=1*u.mJy):
    np.random.seed(42)

    lamb = np.linspace(4000, 10000, 1000)*u.AA
    flux = np.random.randn(1000) * amplitude
    return Spectrum1D(spectral_axis=lamb, flux=flux)


def test_noise_region_unc():
    nspec = pure_noise_spectrum()

    region = SpectralRegion(5000*u.AA, 6000*u.AA)
    unc_spec = noise_region_uncertainty(nspec, region)

    assert np.allclose(unc_spec.uncertainty.array, 1, atol=0.1)


@pytest.mark.xfail
def test_noise_region_unc_multi():
    nspec = pure_noise_spectrum()

    multi_region = SpectralRegion([(5000*u.AA, 6000*u.AA), (900*u.pixel, 980*u.pixel)])

    unc_spec = noise_region_uncertainty(nspec, multi_region)

    assert quantity_allclose(unc_spec.uncertainty, 1*u.mJy, atol=.01)  # this is a guess at the atol... need to finalize with real data


def test_noise_estimate_uncertainty():
    # an alternate test with a more complicated input spectrum

    np.random.seed(42)

    # Create values for spectrum.
    frequencies = np.linspace(1, 100, 10000) * u.um
    g = models.Gaussian1D(amplitude=1*u.Jy, mean=10*u.um, stddev=1*u.um)
    noise = np.random.normal(0., 0.01, frequencies.shape) * u.Jy
    flux = g(frequencies) + noise
    spectrum = Spectrum1D(spectral_axis=frequencies, flux=flux)

    # Assign uncertainties using the default standard deviation
    spectral_region = SpectralRegion(50*u.um, 80*u.um)
    spectrum_with_uncertainty = noise_region_uncertainty(spectrum, spectral_region)

    indices = np.nonzero((frequencies >= 50*u.um) & (frequencies <= 80*u.um))
    expected_uncertainty = np.std(flux[indices])*np.ones(len(frequencies))

    assert quantity_allclose(spectrum_with_uncertainty.uncertainty.array,
                             expected_uncertainty.value)
    assert isinstance(spectrum_with_uncertainty.uncertainty, StdDevUncertainty)

    # Same idea, but now with variance.
    spectrum_with_uncertainty = noise_region_uncertainty(spectrum, spectral_region, np.var)

    indices = np.nonzero((frequencies >= 50*u.um) & (frequencies <= 80*u.um))
    expected_uncertainty = np.var(flux[indices])*np.ones(len(frequencies))

    assert quantity_allclose(spectrum_with_uncertainty.uncertainty.array,
                             expected_uncertainty.value)
    assert isinstance(spectrum_with_uncertainty.uncertainty, VarianceUncertainty)

    # Same idea, but now with inverse variance.
    spectrum_with_uncertainty = noise_region_uncertainty(spectrum, spectral_region,
                                                         lambda x: 1/np.var(x))

    indices = np.nonzero((frequencies >= 50*u.um) & (frequencies <= 80*u.um))
    expected_uncertainty = 1/np.var(flux[indices])*np.ones(len(frequencies))

    assert quantity_allclose(spectrum_with_uncertainty.uncertainty.array,
                             expected_uncertainty.value)
    assert isinstance(spectrum_with_uncertainty.uncertainty, InverseVariance)

    # Now try with something that does not return Std, Var or IVar type of noise estimation
    with pytest.raises(ValueError):
        noise_region_uncertainty(spectrum, spectral_region, lambda x: np.std(x)**3)
