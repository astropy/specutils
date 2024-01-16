import numpy as np
import astropy.units as u
from astropy import constants as const

from astropy.nddata import StdDevUncertainty
from astropy.modeling import models

from ..spectra.spectrum1d import Spectrum1D
from ..analysis import correlation


def test_autocorrelation():
    """
    Test auto correlation
    """
    size = 42

    # Seed np.random so that results are consistent
    np.random.seed(41)

    # Create test spectra
    spec_axis = np.linspace(5000., 5040., num=size) * u.AA
    f1 = np.random.randn(size) * u.Jy
    g1 = models.Gaussian1D(amplitude=30 * u.Jy, mean=5020 * u.AA, stddev=2 * u.AA)
    flux1 = f1 + g1(spec_axis)

    # Observed spectrum must have a rest wavelength value set in.
    spec1 = Spectrum1D(spectral_axis=spec_axis,
                       flux=flux1,
                       uncertainty=StdDevUncertainty(np.random.sample(size), unit='Jy'),
                       velocity_convention='optical',
                       rest_value=5020.*u.AA)

    spec2 = Spectrum1D(spectral_axis=spec_axis,
                       flux=flux1,
                       uncertainty=StdDevUncertainty(np.random.sample(size), unit='Jy'))

    # Get result from correlation
    corr, lag = correlation.template_correlate(spec1, spec2)

    # Check units
    assert corr.unit == u.dimensionless_unscaled
    assert lag.unit == u.km / u.s

    # Check that lags are symmetrical
    midpoint = int(len(lag) / 2)
    np.testing.assert_almost_equal(lag[midpoint+1].value, (-(lag[midpoint-1])).value, 2)

    # Check position of correlation peak.
    maximum = np.argmax(corr)
    assert maximum == midpoint

    # Test that this works without uncertainty
    spec1 = Spectrum1D(spectral_axis=spec_axis,
                       flux=flux1,
                       velocity_convention='optical',
                       rest_value=5020.*u.AA)

    corr, lag = correlation.template_correlate(spec1, spec2)

    assert corr.unit == u.dimensionless_unscaled
    assert lag.unit == u.km / u.s


def test_correlation():
    """
    Test correlation when both observed and template spectra have the same wavelength axis
    """
    size = 4000

    # Seed np.random so that results are consistent
    np.random.seed(51)

    # Create test spectra
    spec_axis = np.linspace(4001., 8000., num=size) * u.AA

    # Two narrow Gaussians are offset from each other so
    # as to generate a correlation peak at a expected lag.
    f1 = np.random.randn(size)*0.5 * u.Jy
    f2 = np.random.randn(size)*0.5 * u.Jy

    expected_lag = 10000. * u.km/u.s
    mean2 = 5000. * u.AA
    stdev2 = 30. * u.AA
    mean1 = (1 + expected_lag / const.c.to('km/s')) * mean2
    stdev1 = (1 + expected_lag / const.c.to('km/s')) * stdev2
    rest_value = mean2

    g1 = models.Gaussian1D(amplitude=30 * u.Jy, mean=mean1, stddev=stdev1)
    g2 = models.Gaussian1D(amplitude=30 * u.Jy, mean=mean2, stddev=stdev2)

    flux1 = f1 + g1(spec_axis)
    flux2 = f2 + g2(spec_axis)

    # Observed spectrum must have a rest wavelength value set in.
    spec1 = Spectrum1D(spectral_axis=spec_axis,
                       flux=flux1,
                       uncertainty=StdDevUncertainty(np.random.sample(size), unit='Jy'),
                       velocity_convention='optical',
                       rest_value=rest_value)

    spec2 = Spectrum1D(spectral_axis=spec_axis,
                       flux=flux2,
                       uncertainty=StdDevUncertainty(np.random.sample(size), unit='Jy'))

    # Get result from correlation
    corr, lag = correlation.template_correlate(spec1, spec2)

    # Check units
    assert corr.unit == u.dimensionless_unscaled
    assert lag.unit == u.km / u.s

    # Check position of correlation peak.
    maximum = np.argmax(corr)
    v_fit = _fit_peak(corr, lag, maximum)
    # checks against 1.5 * 10**(-decimal)
    np.testing.assert_almost_equal(v_fit.value, expected_lag.value, -1)


def _create_arrays(size1, size2):
    spec_axis_1 = np.linspace(6050., 6100., num=size1) * u.AA
    spec_axis_2 = np.linspace(6000., 6100., num=size2) * u.AA

    # Create continuum-subtracted flux arrays with noise.
    f1 = np.random.randn(size1) * u.Jy
    f2 = np.random.randn(size2) * u.Jy

    # Two narrow Gaussians are offset from each other so
    # as to generate a correlation peak at a expected lag.
    expected_lag = 1000. * u.km/u.s
    mean2 = 6050. * u.AA
    stdev2 = 5. * u.AA
    mean1 = (1 + expected_lag / const.c.to('km/s')) * mean2
    stdev1 = (1 + expected_lag / const.c.to('km/s')) * stdev2
    rest_value = mean2

    g1 = models.Gaussian1D(amplitude=50 * u.Jy, mean=mean1, stddev=stdev1)
    g2 = models.Gaussian1D(amplitude=50 * u.Jy, mean=mean2, stddev=stdev2)

    flux1 = f1 + g1(spec_axis_1)
    flux2 = f2 + g2(spec_axis_2)

    return spec_axis_1, spec_axis_2, flux1, flux2, expected_lag, rest_value


def _fit_peak(corr, lag, index_peak):

    # Parabolic fit to maximum

    n = 3  # points to the left and right of correlation maximum

    peak_lags = lag[index_peak - n:index_peak + n + 1].value
    peak_vals = corr[index_peak - n:index_peak + n + 1].value
    p = np.polyfit(peak_lags, peak_vals, deg=2)
    roots = np.roots(p)

    return np.mean(roots) * u.km / u.s  # maximum lies at mid point between roots


def test_correlation_zero_padding():
    """
    Test correlation when observed and template spectra have different sizes
    """
    size1 = 51
    size2 = 101

    # Seed np.random so that results are consistent
    np.random.seed(41)

    # Create test spectra
    spec_axis_1, spec_axis_2, flux1, flux2, expected_lag, rest_value = _create_arrays(size1, size2)

    # Observed spectrum must have a rest wavelength value set in.
    # Uncertainty is arbitrary.
    spec1 = Spectrum1D(spectral_axis=spec_axis_1,
                       flux=flux1,
                       uncertainty=StdDevUncertainty(np.random.sample(size1), unit='Jy'),
                       velocity_convention='optical',
                       rest_value=rest_value)

    spec2 = Spectrum1D(spectral_axis=spec_axis_2,
                       flux=flux2,
                       uncertainty=StdDevUncertainty(np.random.sample(size2), unit='Jy'))

    # Get result from correlation
    corr, lag = correlation.template_correlate(spec1, spec2)

    # Check units
    assert corr.unit == u.dimensionless_unscaled
    assert lag.unit == u.km / u.s

    # Check that lags are symmetrical
    midpoint = int(len(lag) / 2)
    np.testing.assert_almost_equal(lag[midpoint+1].value, (-(lag[midpoint-1])).value, 2)

    # Check position of correlation peak.
    maximum = np.argmax(corr)
    v_fit = _fit_peak(corr, lag, maximum)
    # checks against 1.5 * 10**(-decimal)
    np.testing.assert_almost_equal(v_fit.value, expected_lag.value, -1)


def test_correlation_random_lines():
    """
    Test correlation when observed and template spectra have different sizes, and
    there are non-correlated features in the observed and template spectra.
    """
    size1 = 51
    size2 = 101

    # Seed np.random so that results are consistent
    np.random.seed(41)

    # Create test spectra
    spec_axis_1, spec_axis_2, flux1, flux2, expected_lag, rest_value = _create_arrays(size1, size2)

    # Add random lines to both spectra to simulate non-correlated spectral features.
    #
    # The more lines we add, and the larger amplitude they have, the larger the error in correlation
    # peak position will be. In this specific case, increasing nlines by 1, or increasing the
    # Gaussian amplitudes by a bit, is enough to make the correlation peak move off position.
    #
    # So we can say generically (as expected anyway) that the presence of non-correlated features
    # in observed and template spectra generates an error in the correlation peak position, and
    # that error will be larger as these non-correlated features become more predominant in the
    # data.
    nlines = 20
    for i in range(nlines):
        mean = (spec_axis_1[-1] - spec_axis_1[0]) * np.random.randn(size1) + spec_axis_1[0]
        g1 = models.Gaussian1D(amplitude=10 * u.Jy, mean=mean, stddev=4 * u.AA)
        flux1 += g1(spec_axis_1)
    for i in range(nlines):
        mean = (spec_axis_2[-1] - spec_axis_2[0]) * np.random.randn(size2) + spec_axis_2[0]
        g2 = models.Gaussian1D(amplitude=5 * u.Jy, mean=mean, stddev=4 * u.AA)
        flux2 += g2(spec_axis_2)

    # Observed spectrum must have a rest wavelength value set in.
    # Uncertainty is arbitrary.
    spec1 = Spectrum1D(spectral_axis=spec_axis_1,
                       flux=flux1,
                       uncertainty=StdDevUncertainty(np.random.sample(size1), unit='Jy'),
                       velocity_convention='optical',
                       rest_value=rest_value)

    spec2 = Spectrum1D(spectral_axis=spec_axis_2,
                       flux=flux2,
                       uncertainty=StdDevUncertainty(np.random.sample(size2), unit='Jy'))

    # Get result from correlation
    corr, lag = correlation.template_correlate(spec1, spec2)

    # Check units
    assert corr.unit == u.dimensionless_unscaled
    assert lag.unit == u.km / u.s

    # Check that lags are symmetrical
    midpoint = int(len(lag) / 2)
    np.testing.assert_almost_equal(lag[midpoint+1].value, (-(lag[midpoint-1])).value, 2)

    # Check position of correlation peak.
    maximum = np.argmax(corr)
    v_fit = _fit_peak(corr, lag, maximum)
    # checks against 1.5 * 10**(-decimal)
    np.testing.assert_almost_equal(v_fit.value, expected_lag.value, -1)
