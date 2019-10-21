import astropy.units as u
import numpy as np
from astropy.nddata import StdDevUncertainty
from astropy.modeling import models
import astropy.units as u

from ..spectra.spectrum1d import Spectrum1D
from ..analysis import correlation

SIZE = 40


def test_correlation():
    """
    Test correlation when both observed and template spectra have the same wavelength axis
    """
    # Seed np.random so that results are consistent
    np.random.seed(41)

    # Create test spectra
    spec_axis = np.linspace(0, SIZE, SIZE) * u.AA

    # Two narrow Gaussians are offset from each other so
    # as to generate a correlation peak at a expected lag.
    f1 = np.random.randn(SIZE) * u.Jy
    f2 = np.random.randn(SIZE) * u.Jy
    g1 = models.Gaussian1D(amplitude=30 * u.Jy, mean=20 * u.AA, stddev=2 * u.AA)
    g2 = models.Gaussian1D(amplitude=30 * u.Jy, mean=24 * u.AA, stddev=2 * u.AA)

    flux1 = f1 + g1(spec_axis)
    flux2 = f2 + g2(spec_axis)

    spec1 = Spectrum1D(spectral_axis=spec_axis,
                      flux=flux1,
                      uncertainty=StdDevUncertainty(np.random.sample(SIZE), unit='Jy'))

    spec2 = Spectrum1D(spectral_axis=spec_axis,
                       flux=flux2,
                       uncertainty=StdDevUncertainty(np.random.sample(SIZE), unit='Jy'))

    # Get result from correlation
    corr_result = correlation.template_correlate(spec1, spec2)

    # Check that lag at mid-point is zero and lags are symmetrical
    assert int((corr_result.spectral_axis[int(len(corr_result.spectral_axis)/2)]).value) == 0
    assert corr_result.spectral_axis[0] == -(corr_result.spectral_axis[-1])

    # Check units
    assert corr_result.flux.unit == u.dimensionless_unscaled
    assert corr_result.spectral_axis.unit == spec1.spectral_axis.unit

    # Lag step must be identical with wavelength step.
    lag = corr_result.spectral_axis[1] - corr_result.spectral_axis[0]
    wstep = spec1.spectral_axis[1] - spec1.spectral_axis[0]
    np.testing.assert_almost_equal(lag.value, wstep.value, 0.0001)

    # Check position of correlation peak.
    assert np.argmax(corr_result.flux) == 35


