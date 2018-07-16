from astropy.modeling import models, fitting
import numpy as np

from specutils.tests.spectral_examples import simulated_spectra
from specutils.spectra.spectrum1d import Spectrum1D
from specutils.fitting import fit_lines


def single_peak():
    np.random.seed(0)
    x = np.linspace(0., 10., 200)
    y_single = 3 * np.exp(-0.5 * (x - 6.3)**2 / 0.8**2)
    y_single += np.random.normal(0., 0.2, x.shape)
    return x, y_single


def single_peak_continuum():
    np.random.seed(0)
    x = np.linspace(0., 10., 200)
    y_single = 3 * np.exp(-0.5 * (x - 6.3)**2 / 0.3**2)
    y_single += np.random.normal(0., 0.2, x.shape)

    y_continuum = 3.2 * np.exp(-0.5 * (x - 0.6)**2 / 2.8**2)
    y_single += y_continuum
    return x, y_single


def single_peak_extra():
    x, y_single = single_peak()
    extra = 4 * np.exp(-0.5 * (x + 8.3)**2 / 0.1**2)
    y_single_extra = y_single + extra
    return x, y_single_extra


def double_peak():
    np.random.seed(42)
    g1 = models.Gaussian1D(1, 5, 0.2)
    g2 = models.Gaussian1D(2.5, 5.5, 0.1)
    x = np.linspace(0, 10, 200)
    y_double = g1(x) + g2(x) + np.random.normal(0., 0.2, x.shape)
    return x, y_double

def test_single_peak_fit():
    """
    Single Peak fit.
    """

    x_single, y_single = single_peak()
    s_single = Spectrum1D(flux=y_single, spectral_axis=x_single)

    g_init = models.Gaussian1D(amplitude=3., mean=5.5, stddev=1.)
    g_fit = fit_lines(s_single, g_init)
    y_single_fit = g_fit(x_single)

    # Comparing every 10th value.
    y_single_fit_expected = np.array([3.69669474e-13, 3.57992454e-11, 2.36719426e-09, 1.06879318e-07,
               3.29498310e-06, 6.93605383e-05, 9.96945607e-04, 9.78431032e-03,
               6.55675141e-02, 3.00017760e-01, 9.37356842e-01, 1.99969007e+00,
               2.91286375e+00, 2.89719280e+00, 1.96758892e+00, 9.12412206e-01,
               2.88900005e-01, 6.24602556e-02, 9.22061121e-03, 9.29427266e-04])

    assert np.allclose(y_single_fit[::10], y_single_fit_expected, atol=1e-3)


def test_single_peak_fit_window():
    """
    Single Peak fit with a window specified
    """

    x_single, y_single = single_peak()
    s_single = Spectrum1D(flux=y_single, spectral_axis=x_single)

    g_init = models.Gaussian1D(amplitude=3., mean=5.5, stddev=1.)
    g_fit = fit_lines(s_single, g_init, window=2)
    y_single_fit = g_fit(x_single)

    # Comparing every 10th value.
    y_single_fit_expected = np.array([5.37829623e-13, 4.90533260e-11, 3.07071391e-09, 1.31934207e-07,
           3.89065408e-06, 7.87471756e-05, 1.09394253e-03, 1.04304054e-02,
           6.82582302e-02, 3.06588610e-01, 9.45157456e-01, 1.99985943e+00,
           2.90430306e+00, 2.89488656e+00, 1.98047026e+00, 9.29934277e-01,
           2.99697654e-01, 6.62920712e-02, 1.00643791e-02, 1.04871999e-03])

    assert np.allclose(y_single_fit[::10], y_single_fit_expected, atol=1e-3)


def test_single_peak_fit_tuple_window():
    """
    Single Peak fit with a window specified as a tuple
    """

    x_single, y_single = single_peak()
    s_single = Spectrum1D(flux=y_single, spectral_axis=x_single)

    g_init = models.Gaussian1D(amplitude=3., mean=5.5, stddev=1.)
    g_fit = fit_lines(s_single, g_init, window=(6, 7))
    y_single_fit = g_fit(x_single)

    # Comparing every 10th value.
    y_single_fit_expected = np.array([2.27910247e-16, 6.61153507e-14, 1.19929786e-11, 1.36031069e-09,
           9.64795926e-08, 4.27877165e-06, 1.18655824e-04, 2.05752587e-03,
           2.23093979e-02, 1.51257422e-01, 6.41256583e-01, 1.69993877e+00,
           2.81787042e+00, 2.92075452e+00, 1.89302084e+00, 7.67188406e-01,
           1.94417326e-01, 3.08073399e-02, 3.05252853e-03, 1.89126142e-04])

    assert np.allclose(y_single_fit[::10], y_single_fit_expected, atol=1e-3)


def test_double_peak_fit():
    """
    Double Peak fit.
    """

    x_double, y_double = double_peak()
    s_double = Spectrum1D(flux=y_double, spectral_axis=x_double, rest_value=0)

    g1_init = models.Gaussian1D(amplitude=2., mean=5.8, stddev=0.2)
    g2_init = models.Gaussian1D(amplitude=1., mean=4.5, stddev=0.2)
    g12_fit = fit_lines(s_double, g1_init+g2_init)
    y12_double_fit = g12_fit(x_double)

    # Comparing every 10th value.
    y12_double_fit_expected = np.array([6.41267744e-177, 3.76609436e-143, 5.68693690e-113, 2.20800998e-086,
               2.20424095e-063, 5.65786408e-044, 3.73406063e-028, 6.33644444e-016,
               2.76468278e-007, 3.10156148e-002, 8.96121804e-001, 2.22150879e+000,
               1.25211690e-004, 2.72881561e-016, 7.82165893e-031, 2.53533026e-047,
               2.11308760e-067, 4.52830712e-091, 2.49510977e-118, 3.53491398e-149])

    assert np.allclose(y12_double_fit[::10], y12_double_fit_expected, atol=1e-3)


def test_double_peak_fit_tuple_window():
    """
    Doulbe Peak fit with a window specified as a tuple
    """

    x_double, y_double = double_peak()
    s_double = Spectrum1D(flux=y_double, spectral_axis=x_double, rest_value=0)

    g2_init = models.Gaussian1D(amplitude=1., mean=4.7, stddev=0.2)
    g2_fit = fit_lines(s_double, g2_init, window=(4.3, 5.3))
    y2_double_fit = g2_fit(x_double)

    # Comparing every 10th value.
    y2_double_fit_expected = np.array([6.05129659e-113, 1.32725733e-091, 1.63739289e-072, 1.13616588e-055,
               4.43426518e-041, 9.73402942e-029, 1.20186233e-018, 8.34656662e-011,
               3.26025884e-005, 7.16287502e-002, 8.85143615e-001, 6.15221044e-002,
               2.40513633e-005, 5.28858259e-011, 6.54078626e-019, 4.55000316e-029,
               1.78026443e-041, 3.91785542e-056, 4.84957155e-073, 3.37636518e-092])

    assert np.allclose(y2_double_fit[::10], y2_double_fit_expected, atol=1e-3)


def test_double_peak_fit_window():
    """
    Double Peak fit with a window.
    """

    x_double, y_double = double_peak()
    s_double = Spectrum1D(flux=y_double, spectral_axis=x_double, rest_value=0)

    g2_init = models.Gaussian1D(amplitude=1., mean=4.7, stddev=0.2)
    g2_fit = fit_lines(s_double, g2_init, window=0.3)
    y2_double_fit = g2_fit(x_double)

    # Comparing every 10th value.
    y2_double_fit_expected = np.array([1.03687057e-208, 8.51224694e-169, 3.94108693e-133, 1.02905822e-101,
               1.51536066e-074, 1.25847387e-051, 5.89419538e-033, 1.55688661e-018,
               2.31921800e-008, 1.94840077e-002, 9.23139096e-001, 2.46665339e-003,
               3.71707429e-010, 3.15897717e-021, 1.51406271e-036, 4.09254414e-056,
               6.23871307e-080, 5.36350616e-108, 2.60048700e-140, 7.11070597e-177])

    assert np.allclose(y2_double_fit[::10], y2_double_fit_expected, atol=1e-3)


def test_double_peak_fit_separate_window():
    """
    Double Peak fit with a window.
    """

    x_double, y_double = double_peak()
    s_double = Spectrum1D(flux=y_double, spectral_axis=x_double, rest_value=0)

    gl_init = models.Gaussian1D(amplitude=1., mean=4.8, stddev=0.2)
    gr_init = models.Gaussian1D(amplitude=2., mean=5.3, stddev=0.2)
    gl_fit, gr_fit = fit_lines(s_double, [gl_init, gr_init], window=0.2)
    yl_double_fit = gl_fit(x_double)
    yr_double_fit = gr_fit(x_double)

    # Comparing every 10th value.
    yl_double_fit_expected = np.array([2.94747023e-312, 3.99888064e-252, 1.96020745e-198, 3.47168419e-151,
               2.22153403e-110, 5.13618088e-076, 4.29044393e-048, 1.29490746e-026,
               1.41204917e-011, 5.56333850e-003, 7.91946299e-001, 4.07315149e-005,
               7.56902561e-016, 5.08187488e-033, 1.23277138e-056, 1.08047743e-086,
               3.42155430e-123, 3.91476554e-166, 1.61831351e-215, 2.41709735e-271])

    assert np.allclose(yl_double_fit[::10], yl_double_fit_expected, atol=1e-3)

    # Comparing every 10th value.
    yr_double_fit_expected = np.array([6.44733204e-153, 2.26531185e-125, 1.44273676e-100, 1.66554871e-078,
               3.48528842e-059, 1.32199937e-042, 9.08939374e-029, 1.13279083e-017,
               2.55903177e-009, 1.04788294e-003, 7.77787582e-001, 1.04645458e+000,
               2.55206011e-003, 1.12816484e-008, 9.03993610e-017, 1.31301384e-027,
               3.45688068e-041, 1.64972143e-057, 1.42707890e-076, 2.23767235e-098])

    assert np.allclose(yr_double_fit[::10], yr_double_fit_expected, atol=1e-3)


def test_double_peak_fit_separate_window_tuple_window():
    """
    Double Peak fit with a window.
    """

    x_double, y_double = double_peak()
    s_double = Spectrum1D(flux=y_double, spectral_axis=x_double, rest_value=0)

    g1_init = models.Gaussian1D(amplitude=2., mean=5.3, stddev=0.2)
    g2_init = models.Gaussian1D(amplitude=1., mean=4.9, stddev=0.1)
    g1_fit, g2_fit = fit_lines(s_double, [g1_init, g2_init], window=[(5.3, 5.8), (4.6, 5.3)])
    y1_double_fit = g1_fit(x_double)
    y2_double_fit = g2_fit(x_double)

    # Comparing every 10th value.
    y1_double_fit_expected = np.array([0.00000000e+000, 0.00000000e+000, 4.33665720e-275, 6.16153689e-217,
               1.08340069e-165, 2.35751933e-121, 6.34873955e-084, 2.11585478e-053,
               8.72670508e-030, 4.45431410e-013, 2.81369795e-003, 2.19958225e+000,
               2.12798732e-004, 2.54779195e-015, 3.77506978e-033, 6.92232688e-058,
               1.57088835e-089, 4.41168638e-128, 1.53331026e-173, 6.59510245e-226])

    assert np.allclose(y1_double_fit[::10], y1_double_fit_expected, atol=1e-3)

    # Comparing every 10th value.
    y2_double_fit_expected = np.array([0.00000000e+000, 0.00000000e+000, 0.00000000e+000, 1.22905611e-250,
           4.39298703e-182, 1.69508976e-124, 7.06106233e-078, 3.17535260e-042,
           1.54155292e-017, 8.07922093e-004, 4.57114606e-001, 2.79206342e-009,
           1.84106884e-028, 1.31056734e-058, 1.00714783e-099, 8.35548514e-152,
           7.48332638e-215, 7.23539540e-289, 0.00000000e+000, 0.00000000e+000])

    assert np.allclose(y2_double_fit[::10], y2_double_fit_expected, atol=1e-3)


def test_double_peak_fit_with_exclusion():
    """
    Double Peak fit with a window.
    """

    x_double, y_double = double_peak()
    s_double = Spectrum1D(flux=y_double, spectral_axis=x_double, rest_value=0)

    g1_init = models.Gaussian1D(amplitude=1., mean=4.9, stddev=0.2)
    g1_fit = fit_lines(s_double, g1_init, exclude_regions=[(5.2, 5.8)])
    y1_double_fit = g1_fit(x_double)

    # Comparing every 10th value.
    y1_double_fit_expected = np.array([1.93936382e-91, 3.15688868e-74, 7.94634186e-59, 3.09302350e-45,
               1.86168925e-33, 1.73276298e-23, 2.49390054e-15, 5.55043594e-09,
               1.91021855e-04, 1.01659277e-01, 8.36602151e-01, 1.06463057e-01,
               2.09501379e-04, 6.37503889e-09, 2.99976192e-15, 2.18272355e-23,
               2.45594524e-33, 4.27313508e-45, 1.14969479e-58, 4.78328766e-74])

    assert np.allclose(y1_double_fit[::10], y1_double_fit_expected, atol=1e-3)
