import os
import numpy as np
import astropy.units as u

from astropy.nddata import StdDevUncertainty
from astropy.modeling import models

from ..spectra.spectrum1d import Spectrum1D
from ..analysis import correlation

SIZE = 41 # force symmetry around peak.


def test_correlation():
    """
    Test correlation when both observed and template spectra have the same wavelength axis
    """
    # Seed np.random so that results are consistent
    np.random.seed(41)

    # Create test spectra
    spec_axis = np.linspace(5000., 5040., num=SIZE) * u.AA

    # Two narrow Gaussians are offset from each other so
    # as to generate a correlation peak at a expected lag.
    f1 = np.random.randn(SIZE) * u.Jy
    f2 = np.random.randn(SIZE) * u.Jy

    g1 = models.Gaussian1D(amplitude=30 * u.Jy, mean=5025 * u.AA, stddev=2 * u.AA)
    g2 = models.Gaussian1D(amplitude=30 * u.Jy, mean=5020 * u.AA, stddev=2 * u.AA)

    flux1 = f1 + g1(spec_axis)
    flux2 = f2 + g2(spec_axis)

    # Observed spectrum must have a rest wavelength value set in.
    spec1 = Spectrum1D(spectral_axis=spec_axis,
                      flux=flux1,
                      uncertainty=StdDevUncertainty(np.random.sample(SIZE), unit='Jy'),
                      velocity_convention='optical',
                      rest_value=spec_axis[int(SIZE/2)])

    spec2 = Spectrum1D(spectral_axis=spec_axis,
                       flux=flux2,
                       uncertainty=StdDevUncertainty(np.random.sample(SIZE), unit='Jy'))

    # Get result from correlation
    corr, lag = correlation.template_correlate(spec1, spec2)

    # Check units
    assert corr.unit == u.dimensionless_unscaled
    assert lag.unit == u.km / u.s

    # Check position of correlation peak.
    maximum = np.argmax(corr)
    assert maximum == 45
    np.testing.assert_almost_equal(lag[maximum].value, 298.6, 1)

    # Check that lags are symmetrical
    midpoint = int(len(lag) / 2)
    np.testing.assert_almost_equal(lag[midpoint+2].value, (-(lag[midpoint-2])).value, 2)
