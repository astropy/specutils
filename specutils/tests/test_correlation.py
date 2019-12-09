import os
import numpy as np
import astropy.units as u
from astropy import constants as const

from astropy.nddata import StdDevUncertainty
from astropy.modeling import models

from ..spectra.spectrum1d import Spectrum1D
from ..analysis import correlation
from ..manipulation.resample import FluxConservingResampler

SIZE = 1001

DATA_DIR = os.path.join(os.path.dirname(__file__), '../../')


def test_correlation():
    """
    Test correlation when both observed and template spectra have the same wavelength axis
    """
    # Seed np.random so that results are consistent
    np.random.seed(41)

    # Create test spectra
    spec_axis = np.linspace(4500., 6500., num=SIZE) * u.AA

    # Two narrow Gaussians are offset from each other so
    # as to generate a correlation peak at a expected lag.
    # f1 = np.random.randn(SIZE)*0.0001 * u.Jy
    # f2 = np.random.randn(SIZE)*0.0001 * u.Jy
    f1 = (np.ones(SIZE)-0.5) * u.Jy
    f2 = (np.ones(SIZE)-0.5) * u.Jy

    rest_value = 6000. * u.AA

    mean1 = 5035. * u.AA
    mean2 = 5015. * u.AA
    expected_lag = (mean1 - mean2) / rest_value * const.c.to('km/s')

    g1 = models.Gaussian1D(amplitude=30 * u.Jy, mean=mean1, stddev=10. * u.AA)
    g2 = models.Gaussian1D(amplitude=30 * u.Jy, mean=mean2, stddev=10. * u.AA)

    flux1 = f1 + g1(spec_axis)
    flux2 = f2 + g2(spec_axis)

    # Observed spectrum must have a rest wavelength value set in.
    spec1 = Spectrum1D(spectral_axis=spec_axis,
                      flux=flux1,
                      uncertainty=StdDevUncertainty(np.random.sample(SIZE), unit='Jy'),
                      velocity_convention='optical',
                      rest_value=rest_value)

    spec2 = Spectrum1D(spectral_axis=spec_axis,
                       flux=flux2,
                       uncertainty=StdDevUncertainty(np.random.sample(SIZE), unit='Jy'))

    # Get result from correlation
    corr, lag = correlation.template_correlate(spec1, spec2)

    # Check units
    assert corr.unit == u.dimensionless_unscaled
    assert lag.unit == u.km / u.s

    # Check position of correlation peak.
    corr_peak = np.where(corr == np.amax(corr))[0][0]
    np.testing.assert_almost_equal(lag[corr_peak].value, expected_lag.value, 1)
