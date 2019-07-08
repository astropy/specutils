import numpy as np

import astropy.units as u

from ..spectra.spectrum1d import Spectrum1D
from ..fitting.continuum import fit_generic_continuum


def single_peak_continuum():
    np.random.seed(0)
    x = np.linspace(0., 10., 200)
    y_single = 3 * np.exp(-0.5 * (x - 6.3)**2 / 0.1**2)
    y_single += np.random.normal(0., 0.2, x.shape)

    y_continuum = 3.2 * np.exp(-0.5 * (x - 5.6)**2 / 4.8**2)
    y_single += y_continuum
    return x, y_single


def test_continuum_fit():
    """
    This test fits the the first simulated spectrum from the fixture.  The
    initial guesses are manually set here with bounds that essentially make
    sense as the functionality of the test is to make sure the fit works and
    we get a reasonable answer out **given** good initial guesses.
    """

    x_single_continuum, y_single_continuum = single_peak_continuum()
    s_single_continuum = Spectrum1D(flux=y_single_continuum*u.Jy, spectral_axis=x_single_continuum*u.um)
    g1_fit = fit_generic_continuum(s_single_continuum)

    y_continuum_fitted = g1_fit(s_single_continuum.spectral_axis)

    y_continuum_fitted_expected = np.array([1.71364056, 1.87755574, 2.05310622, 2.23545755, 2.41977527,
                                            2.60122493, 2.77497207, 2.93618225, 3.080021, 3.20165388,
                                            3.29624643, 3.3589642, 3.38497273, 3.36943758, 3.30752428,
                                            3.19439839, 3.02522545, 2.79517101, 2.49940062, 2.13307982])

    assert np.allclose(y_continuum_fitted.value[::10], y_continuum_fitted_expected, atol=1e-5)


def test_continuum_calculation():
    """
    This test fits the the first simulated spectrum from the fixture.  The
    initial guesses are manually set here with bounds that essentially make
    sense as the functionality of the test is to make sure the fit works and
    we get a reasonable answer out **given** good initial guesses.
    """

    x_single_continuum, y_single_continuum = single_peak_continuum()
    spectrum = Spectrum1D(flux=y_single_continuum*u.Jy, spectral_axis=x_single_continuum*u.um)
    g1_fit = fit_generic_continuum(spectrum)

    spectrum_normalized = spectrum / g1_fit(spectrum.spectral_axis)

    y_continuum_fitted_expected = np.array([1.15139925, 0.98509363, 0.73700614, 1.00911864, 0.913129,
                                            0.93145533, 0.94904202, 1.04162879, 0.90851397, 0.9494352,
                                            1.07812394, 1.06376489, 0.98705237, 0.94569623, 0.83502377,
                                            0.91909416, 0.89662208, 1.01458511, 0.96124191, 0.94847744])

    assert np.allclose(spectrum_normalized.flux.value[::10], y_continuum_fitted_expected, atol=1e-5)
