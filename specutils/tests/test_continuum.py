import numpy as np

import astropy.units as u

from specutils.spectra.spectrum1d import Spectrum1D
from specutils.fitting.continuum import fit_continuum_generic


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
    g1_fit = fit_continuum_generic(s_single_continuum)

    y_continuum_fitted = g1_fit(x_single_continuum*u.um)

    y_continuum_fitted_expected = np.array([1.71414049, 1.87778562, 2.05313605, 2.23534949, 2.41958364,
               2.60099619, 2.77474484, 2.93598729, 3.07988123, 3.20158436,
               3.29625438, 3.35904898, 3.38512587, 3.36964273, 3.30775726,
               3.19462717, 3.02541014, 2.79526388, 2.49934609, 2.13281445])

    assert np.allclose(y_continuum_fitted.value[::10], y_continuum_fitted_expected, atol=1e-5)
