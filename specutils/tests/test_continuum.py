import numpy as np
import pytest

import astropy.units as u
from astropy.modeling.polynomial import Chebyshev1D
from astropy.utils.exceptions import AstropyUserWarning

from ..spectra.spectrum1d import Spectrum1D
from ..spectra import SpectralRegion
from ..fitting.continuum import fit_generic_continuum, fit_continuum
from ..manipulation.smoothing import median_smooth


def single_peak_continuum(noise=0.2, constant_continuum=False):
    np.random.seed(0)
    x = np.linspace(0., 10., 200)
    y_single = 3 * np.exp(-0.5 * (x - 6.3)**2 / 0.1**2)
    y_single += np.random.normal(0., noise, x.shape)
    if not constant_continuum:
        y_continuum = 3.2 * np.exp(-0.5 * (x - 5.6)**2 / 4.8**2)
    else:
        y_continuum = 3.2
    y_single += y_continuum
    return x, y_single


def test_continuum_fit():
    """
    This test fits the first simulated spectrum from the fixture.  The
    initial guesses are manually set here with bounds that essentially make
    sense as the functionality of the test is to make sure the fit works and
    we get a reasonable answer out **given** good initial guesses.
    """

    x_single_continuum, y_single_continuum = single_peak_continuum()
    s_single_continuum = Spectrum1D(flux=y_single_continuum*u.Jy, spectral_axis=x_single_continuum*u.um)
    with pytest.warns(AstropyUserWarning, match='Model is linear in parameters'):
        g1_fit = fit_generic_continuum(s_single_continuum)

    y_continuum_fitted = g1_fit(s_single_continuum.spectral_axis)

    y_continuum_fitted_expected = np.array([1.71364056, 1.87755574, 2.05310622, 2.23545755, 2.41977527,
                                            2.60122493, 2.77497207, 2.93618225, 3.080021, 3.20165388,
                                            3.29624643, 3.3589642, 3.38497273, 3.36943758, 3.30752428,
                                            3.19439839, 3.02522545, 2.79517101, 2.49940062, 2.13307982])

    assert np.allclose(y_continuum_fitted.value[::10], y_continuum_fitted_expected, atol=1e-5)


def test_continuum_calculation():
    """
    This test fits the first simulated spectrum from the fixture.  The
    initial guesses are manually set here with bounds that essentially make
    sense as the functionality of the test is to make sure the fit works and
    we get a reasonable answer out **given** good initial guesses.
    """

    x_single_continuum, y_single_continuum = single_peak_continuum()
    spectrum = Spectrum1D(flux=y_single_continuum*u.Jy, spectral_axis=x_single_continuum*u.um)
    with pytest.warns(AstropyUserWarning, match='Model is linear in parameters'):
        g1_fit = fit_generic_continuum(spectrum)

    spectrum_normalized = spectrum / g1_fit(spectrum.spectral_axis)

    y_continuum_fitted_expected = np.array([1.15139925, 0.98509363, 0.73700614, 1.00911864, 0.913129,
                                            0.93145533, 0.94904202, 1.04162879, 0.90851397, 0.9494352,
                                            1.07812394, 1.06376489, 0.98705237, 0.94569623, 0.83502377,
                                            0.91909416, 0.89662208, 1.01458511, 0.96124191, 0.94847744])

    assert np.allclose(spectrum_normalized.flux.value[::10], y_continuum_fitted_expected, atol=1e-5)


def test_continuum_full_window():
    """
    This test fits the first simulated spectrum from the fixture, but
    with the fit_continuum function instead of fit_generic_continuum. Uses
    a window to select the entire spectrum and checks that it recovers the
    original, non-windowed fit.
    """

    x_single_continuum, y_single_continuum = single_peak_continuum()
    spectrum = Spectrum1D(flux=y_single_continuum*u.Jy, spectral_axis=x_single_continuum*u.um)

    # Smooth in the same way fit_generic_continuum does.
    spectrum_smoothed = median_smooth(spectrum, 3)

    # Check that a full width window recovers the original, non-windowed fit.
    with pytest.warns(AstropyUserWarning, match='Model is linear in parameters'):
        g1_fit = fit_continuum(spectrum_smoothed, window=(0.*u.um, 10.*u.um))
    with pytest.warns(AstropyUserWarning, match='Model is linear in parameters'):
        g1_fit_orig = fit_continuum(spectrum_smoothed)

    sp_normalized = spectrum / g1_fit(spectrum.spectral_axis)
    sp_normalized_orig = spectrum / g1_fit_orig(spectrum.spectral_axis)

    assert np.allclose(sp_normalized.flux.value, sp_normalized_orig.flux.value, atol=1e-5)


def test_continuum_spectral_region():
    """
    As before, but with a SpectralRegion as window.
    """

    x_single_continuum, y_single_continuum = single_peak_continuum()
    spectrum = Spectrum1D(flux=y_single_continuum*u.Jy, spectral_axis=x_single_continuum*u.um)

    # Smooth in the same way fit_generic_continuum does.
    spectrum_smoothed = median_smooth(spectrum, 3)

    # Check that a full width window recovers the original, non-windowed fit.
    region = SpectralRegion(0.*u.um, 10.*u.um)
    with pytest.warns(AstropyUserWarning, match='Model is linear in parameters'):
        g1_fit = fit_continuum(spectrum_smoothed, window=region)

    spectrum_normalized = spectrum / g1_fit(spectrum.spectral_axis)

    y_continuum_fitted_expected = np.array([1.15139925, 0.98509363, 0.73700614, 1.00911864, 0.913129,
                                            0.93145533, 0.94904202, 1.04162879, 0.90851397, 0.9494352,
                                            1.07812394, 1.06376489, 0.98705237, 0.94569623, 0.83502377,
                                            0.91909416, 0.89662208, 1.01458511, 0.96124191, 0.94847744])

    assert np.allclose(spectrum_normalized.flux.value[::10], y_continuum_fitted_expected, atol=1e-5)


def test_continuum_window_no_noise():
    """
    This test repeats the setup from the previous tests, but uses windows to
    select specific regions to fit. It uses a synthetic spectrum with no
    noise component. That way, numerical effects coming from the model fit
    itself can be more easily spotted.
    """

    x_single_continuum, y_single_continuum = single_peak_continuum(noise=0.)
    spectrum = Spectrum1D(flux=y_single_continuum*u.Jy, spectral_axis=x_single_continuum*u.um)

    # Smooth in the same way fit_generic_continuum does.
    spectrum_smoothed = median_smooth(spectrum, 3)

    # Window selects the first half of the spectrum.
    with pytest.warns(AstropyUserWarning, match='Model is linear in parameters'):
        g1_fit = fit_continuum(spectrum_smoothed, window=(0.*u.um, 5.*u.um))

    spectrum_normalized = spectrum / g1_fit(spectrum.spectral_axis)

    y_continuum_fitted_expected = np.ones(shape=(spectrum_normalized.spectral_axis.shape))

    # Check fit over the first half. The agreement is not so good as before,
    # probably due to the narrow component at 6.5 um.
    assert np.allclose(spectrum_normalized.flux.value[0:100], y_continuum_fitted_expected[0:100],
                       atol=5.5e-4)

    # Window selects the red end of the spectrum.
    with pytest.warns(AstropyUserWarning, match='Model is linear in parameters'):
        g1_fit = fit_continuum(spectrum_smoothed, window=(8.*u.um, 10.*u.um))

    spectrum_normalized = spectrum / g1_fit(spectrum.spectral_axis)

    # Check fit over the red end.
    assert np.allclose(spectrum_normalized.flux.value[160:], y_continuum_fitted_expected[160:],
                       atol=1.e-5)


def test_double_continuum_window():
    """
    Fit to no-noise spectrum comprised of a exponential continuum plus an emission Gaussian
    """
    x_single_continuum, y_single_continuum = single_peak_continuum(noise=0.,constant_continuum=False)
    spectrum = Spectrum1D(flux=y_single_continuum*u.Jy, spectral_axis=x_single_continuum*u.um)

    # Smooth in the same way fit_generic_continuum does.
    spectrum_smoothed = median_smooth(spectrum, 3)

    # Test spectrum is comprised of an exponential continuum that increases sharply at the
    # long wavelength end, plus a large-amplitude narrow Gaussian. We define the fitting
    # window as a SpectralRegion with sub-regions at the blue and red halves of the spectrum,
    # avoiding the Gaussian in between. The polynomial degree is high to accomodate the large
    # amplitude range of the expoinential continuum.
    region = SpectralRegion([(0. * u.um, 5. * u.um), (8. * u.um, 10. * u.um)])
    with pytest.warns(AstropyUserWarning, match='Model is linear in parameters'):
        g1_fit = fit_continuum(spectrum_smoothed, model=Chebyshev1D(7), window=region)

    spectrum_normalized = spectrum / g1_fit(spectrum.spectral_axis)

    y_continuum_fitted_expected = np.ones(shape=(spectrum_normalized.spectral_axis.shape))

    # Check fit over the windowed regions. Note that we fit a noiseless spectrum so we can
    # actually measure the degree of mismatch between model and data (polynomial X exponential).
    assert np.allclose(spectrum_normalized.flux.value[0:90], y_continuum_fitted_expected[0:90],
                       atol=1.e-4)
    assert np.allclose(spectrum_normalized.flux.value[170:], y_continuum_fitted_expected[170:],
                       atol=1.e-4)


def test_double_continuum_window_alternate():
    """
    Fit to no-noise spectrum comprised of a exponential continuum plus an emission Gaussian

    This is the same test as above, but with a different form for the double window
    specification.

    """
    x_single_continuum, y_single_continuum = single_peak_continuum(noise=0.,constant_continuum=False)
    spectrum = Spectrum1D(flux=y_single_continuum*u.Jy, spectral_axis=x_single_continuum*u.um)

    # Smooth in the same way fit_generic_continuum does.
    spectrum_smoothed = median_smooth(spectrum, 3)

    # Test spectrum is comprised of an exponential continuum that increases sharply at the
    # long wavelength end, plus a large-amplitude narrow Gaussian. We define the fitting
    # window as a SpectralRegion with sub-regions at the blue and red halves of the spectrum,
    # avoiding the Gaussian in between. The polynomial degree is high to accomodate the large
    # amplitude range of the expoinential continuum.
    region = [(0.*u.um, 5.*u.um), (8.* u.um, 10.* u.um)]
    with pytest.warns(AstropyUserWarning, match='Model is linear in parameters'):
        g1_fit = fit_continuum(spectrum_smoothed, model=Chebyshev1D(7), window=region)

    spectrum_normalized = spectrum / g1_fit(spectrum.spectral_axis)

    y_continuum_fitted_expected = np.ones(shape=(spectrum_normalized.spectral_axis.shape))

    # Check fit over the windowed regions. Note that we fit a noiseless spectrum so we can
    # actually measure the degree of mismatch between model and data (polynomial X exponential).
    assert np.allclose(spectrum_normalized.flux.value[0:90], y_continuum_fitted_expected[0:90],
                       atol=1.e-4)
    assert np.allclose(spectrum_normalized.flux.value[170:], y_continuum_fitted_expected[170:],
                       atol=1.e-4)
