import numpy as np

import astropy.units as u
from astropy.modeling import models

from ..spectra import Spectrum1D, SpectralRegion
from ..fitting import fit_lines


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
    g1 = models.Gaussian1D(1, 4.6, 0.2)
    g2 = models.Gaussian1D(2.5, 5.5, 0.1)
    x = np.linspace(0, 10, 200)
    y_double = g1(x) + g2(x) + np.random.normal(0., 0.2, x.shape)
    return x, y_double


def test_single_peak_fit():
    """
    Single Peak fit.
    """

    # Create the spectrum
    x_single, y_single = single_peak()
    s_single = Spectrum1D(flux=y_single*u.Jy, spectral_axis=x_single*u.um)

    # Fit the spectrum
    g_init = models.Gaussian1D(amplitude=3.*u.Jy, mean=6.1*u.um, stddev=1.*u.um)
    g_fit = fit_lines(s_single, g_init)
    y_single_fit = g_fit(x_single*u.um)

    # Comparing every 10th value.
    y_single_fit_expected = np.array([3.69669474e-13, 3.57992454e-11, 2.36719426e-09, 1.06879318e-07,
               3.29498310e-06, 6.93605383e-05, 9.96945607e-04, 9.78431032e-03,
               6.55675141e-02, 3.00017760e-01, 9.37356842e-01, 1.99969007e+00,
               2.91286375e+00, 2.89719280e+00, 1.96758892e+00, 9.12412206e-01,
               2.88900005e-01, 6.24602556e-02, 9.22061121e-03, 9.29427266e-04]) * u.Jy

    assert np.allclose(y_single_fit.value[::10], y_single_fit_expected.value, atol=1e-5)


def test_single_peak_fit_window():
    """
    Single Peak fit with a window specified
    """

    # Create the sepctrum
    x_single, y_single = single_peak()
    s_single = Spectrum1D(flux=y_single*u.Jy, spectral_axis=x_single*u.um)

    # Fit the spectrum
    g_init = models.Gaussian1D(amplitude=3.*u.Jy, mean=5.5*u.um, stddev=1.*u.um)
    g_fit = fit_lines(s_single, g_init, window=2*u.um)
    y_single_fit = g_fit(x_single*u.um)

    # Comparing every 10th value.
    y_single_fit_expected = np.array([3.69669474e-13, 3.57992454e-11, 2.36719426e-09, 1.06879318e-07,
                                      3.29498310e-06, 6.93605383e-05, 9.96945607e-04, 9.78431032e-03,
                                      6.55675141e-02, 3.00017760e-01, 9.37356842e-01, 1.99969007e+00,
                                      2.91286375e+00, 2.89719280e+00, 1.96758892e+00, 9.12412206e-01,
                                      2.88900005e-01, 6.24602556e-02, 9.22061121e-03, 9.29427266e-04]) * u.Jy

    assert np.allclose(y_single_fit.value[::10], y_single_fit_expected.value, atol=1e-5)


def test_single_peak_fit_tuple_window():
    """
    Single Peak fit with a window specified as a tuple
    """

    # Create the spectrum to fit
    x_single, y_single = single_peak()
    s_single = Spectrum1D(flux=y_single*u.Jy, spectral_axis=x_single*u.um)

    # Fit the spectrum
    g_init = models.Gaussian1D(amplitude=3.*u.Jy, mean=5.5*u.um, stddev=1.*u.um)
    g_fit = fit_lines(s_single, g_init, window=(6*u.um, 7*u.um))
    y_single_fit = g_fit(x_single*u.um)

    # Comparing every 10th value.
    y_single_fit_expected = np.array([2.29674788e-16, 6.65518998e-14, 1.20595958e-11, 1.36656472e-09,
                                      9.68395624e-08, 4.29141576e-06, 1.18925100e-04, 2.06096976e-03,
                                      2.23354585e-02, 1.51371211e-01, 6.41529836e-01, 1.70026100e+00,
                                      2.81799025e+00, 2.92071068e+00, 1.89305291e+00, 7.67294570e-01,
                                      1.94485245e-01, 3.08273612e-02, 3.05570344e-03, 1.89413625e-04])*u.Jy

    assert np.allclose(y_single_fit.value[::10], y_single_fit_expected.value, atol=1e-5)


def test_double_peak_fit():
    """
    Double Peak fit.
    """

    # Create the spectrum to fit
    x_double, y_double = double_peak()
    s_double = Spectrum1D(flux=y_double*u.Jy, spectral_axis=x_double*u.um)

    # Fit the spectrum
    g1_init = models.Gaussian1D(amplitude=2.3*u.Jy, mean=5.6*u.um, stddev=0.1*u.um)
    g2_init = models.Gaussian1D(amplitude=1.*u.Jy, mean=4.4*u.um, stddev=0.1*u.um)
    g12_fit = fit_lines(s_double, g1_init+g2_init)
    y12_double_fit = g12_fit(x_double*u.um)

    # Comparing every 10th value.
    y12_double_fit_expected = np.array([2.86790780e-130, 2.12984643e-103, 1.20060032e-079, 5.13707226e-059,
                                        1.66839912e-041, 4.11292970e-027, 7.69608184e-016, 1.09308800e-007,
                                        1.17844042e-002, 9.64333366e-001, 6.04322205e-002, 2.22653307e+000,
                                        5.51964567e-005, 8.13581859e-018, 6.37320251e-038, 8.85834856e-055,
                                        1.05230522e-074, 9.48850399e-098, 6.49412764e-124, 3.37373489e-153])

    assert np.allclose(y12_double_fit.value[::10], y12_double_fit_expected, atol=1e-5)


def test_double_peak_fit_tuple_window():
    """
    Doulbe Peak fit with a window specified as a tuple
    """

    # Create the spectrum to fit
    x_double, y_double = double_peak()
    s_double = Spectrum1D(flux=y_double*u.Jy, spectral_axis=x_double*u.um, rest_value=0*u.um)

    # Fit the spectrum.
    g2_init = models.Gaussian1D(amplitude=1.*u.Jy, mean=4.7*u.um, stddev=0.2*u.um)
    g2_fit = fit_lines(s_double, g2_init, window=(4.3*u.um, 5.3*u.um))
    y2_double_fit = g2_fit(x_double*u.um)

    # Comparing every 10th value.
    y2_double_fit_expected = np.array([2.82386634e-116, 2.84746284e-092, 4.63895634e-071, 1.22104254e-052,
                                       5.19265653e-037, 3.56776869e-024, 3.96051875e-014, 7.10322789e-007,
                                       2.05829545e-002, 9.63624806e-001, 7.28880815e-002, 8.90744929e-006,
                                       1.75872724e-012, 5.61037526e-022, 2.89156942e-034, 2.40781783e-049,
                                       3.23938019e-067, 7.04122962e-088, 2.47276807e-111, 1.40302869e-137])

    assert np.allclose(y2_double_fit.value[::10], y2_double_fit_expected, atol=1e-5)


def test_double_peak_fit_window():
    """
    Double Peak fit with a window.
    """

    # Create the specturm to fit
    x_double, y_double = double_peak()
    s_double = Spectrum1D(flux=y_double*u.Jy, spectral_axis=x_double*u.um, rest_value=0*u.um)

    # Fit the spectrum
    g2_init = models.Gaussian1D(amplitude=1.*u.Jy, mean=4.7*u.um, stddev=0.2*u.um)
    g2_fit = fit_lines(s_double, g2_init, window=0.3*u.um)
    y2_double_fit = g2_fit(x_double*u.um)

    # Comparing every 10th value.
    y2_double_fit_expected = np.array([1.66363393e-128, 5.28910721e-102, 1.40949521e-078, 3.14848385e-058,
                                       5.89516506e-041, 9.25224449e-027, 1.21718016e-015, 1.34220626e-007,
                                       1.24062432e-002, 9.61209273e-001, 6.24240938e-002, 3.39815491e-006,
                                       1.55056770e-013, 5.93054936e-024, 1.90132233e-037, 5.10943886e-054,
                                       1.15092572e-073, 2.17309153e-096, 3.43926290e-122, 4.56256813e-151])

    assert np.allclose(y2_double_fit.value[::10], y2_double_fit_expected, atol=1e-5)


def test_double_peak_fit_separate_window():
    """
    Double Peak fit with a window.
    """

    # Create the spectrum to fit
    x_double, y_double = double_peak()
    s_double = Spectrum1D(flux=y_double*u.Jy, spectral_axis=x_double*u.um, rest_value=0*u.um)

    # Fit the spectrum
    gl_init = models.Gaussian1D(amplitude=1.*u.Jy, mean=4.8*u.um, stddev=0.2*u.um)
    gr_init = models.Gaussian1D(amplitude=2.*u.Jy, mean=5.3*u.um, stddev=0.2*u.um)
    gl_fit, gr_fit = fit_lines(s_double, [gl_init, gr_init], window=0.2*u.um)
    yl_double_fit = gl_fit(x_double*u.um)
    yr_double_fit = gr_fit(x_double*u.um)

    # Comparing every 10th value.
    yl_double_fit_expected = np.array([3.40725147e-18, 5.05500395e-15, 3.59471319e-12, 1.22527176e-09,
                                       2.00182467e-07, 1.56763547e-05, 5.88422893e-04, 1.05866724e-02,
                                       9.12966452e-02, 3.77377148e-01, 7.47690410e-01, 7.10057397e-01,
                                       3.23214276e-01, 7.05201207e-02, 7.37498248e-03, 3.69687164e-04,
                                       8.88245844e-06, 1.02295712e-07, 5.64686114e-10, 1.49410879e-12])

    assert np.allclose(yl_double_fit.value[::10], yl_double_fit_expected, atol=1e-5)

    # Comparing every 10th value.
    yr_double_fit_expected = np.array([0.00000000e+000, 0.00000000e+000, 0.00000000e+000, 3.04416285e-259,
                                       3.85323221e-198, 2.98888589e-145, 1.42075875e-100, 4.13864520e-064,
                                       7.38793226e-036, 8.08191847e-016, 5.41792361e-004, 2.22575901e+000,
                                       5.60338234e-005, 8.64468603e-018, 8.17287853e-039, 4.73508430e-068,
                                       1.68115300e-105, 3.65774659e-151, 4.87693358e-205, 3.98480359e-267])

    assert np.allclose(yr_double_fit.value[::10], yr_double_fit_expected, atol=1e-5)


def test_double_peak_fit_separate_window_tuple_window():
    """
    Double Peak fit with a window.
    """

    x_double, y_double = double_peak()
    s_double = Spectrum1D(flux=y_double*u.Jy, spectral_axis=x_double*u.um, rest_value=0*u.um)

    g1_init = models.Gaussian1D(amplitude=2.*u.Jy, mean=5.3*u.um, stddev=0.2*u.um)
    g2_init = models.Gaussian1D(amplitude=1.*u.Jy, mean=4.9*u.um, stddev=0.1*u.um)
    g1_fit, g2_fit = fit_lines(s_double, [g1_init, g2_init], window=[(5.3*u.um, 5.8*u.um), (4.6*u.um, 5.3*u.um)])
    y1_double_fit = g1_fit(x_double*u.um)
    y2_double_fit = g2_fit(x_double*u.um)

    # Comparing every 10th value.
    y1_double_fit_expected = np.array([0.00000000e+000, 0.00000000e+000, 5.61595149e-307, 3.38362505e-242,
                                       4.27358433e-185, 1.13149721e-135, 6.28008984e-094, 7.30683649e-060,
                                       1.78214929e-033, 9.11192086e-015, 9.76623021e-004, 2.19429562e+000,
                                       1.03350951e-004, 1.02043415e-016, 2.11206194e-036, 9.16388177e-064,
                                       8.33495900e-099, 1.58920023e-141, 6.35191874e-192, 5.32209240e-250])


    assert np.allclose(y1_double_fit.value[::10], y1_double_fit_expected, atol=1e-5)

    # Comparing every 10th value.
    y2_double_fit_expected = np.array([2.52990802e-158, 5.15446435e-126, 2.07577138e-097, 1.65231432e-072,
                                       2.59969849e-051, 8.08482210e-034, 4.96975664e-020, 6.03833143e-010,
                                       1.45016006e-003, 6.88386116e-001, 6.45900222e-002, 1.19788723e-006,
                                       4.39120391e-015, 3.18176751e-027, 4.55691000e-043, 1.28999976e-062,
                                       7.21815119e-086, 7.98324559e-113, 1.74521997e-143, 7.54115780e-178])

    assert np.allclose(y2_double_fit.value[::10], y2_double_fit_expected, atol=1e-3)


def test_double_peak_fit_with_exclusion():
    """
    Double Peak fit with a window.
    """

    x_double, y_double = double_peak()
    s_double = Spectrum1D(flux=y_double*u.Jy, spectral_axis=x_double*u.um, rest_value=0*u.um)

    g1_init = models.Gaussian1D(amplitude=1.*u.Jy, mean=4.9*u.um, stddev=0.2*u.um)
    g1_fit = fit_lines(s_double, g1_init, exclude_regions=[SpectralRegion(5.2*u.um, 5.8*u.um)])
    y1_double_fit = g1_fit(x_double*u.um)

    # Comparing every 10th value.
    y1_double_fit_expected = np.array([4.64465938e-130, 3.11793334e-103, 1.60765691e-079, 6.36698036e-059,
                                       1.93681098e-041, 4.52537486e-027, 8.12148549e-016, 1.11951515e-007,
                                       1.18532671e-002, 9.63961653e-001, 6.02136613e-002, 2.88897581e-006,
                                       1.06464879e-013, 3.01357787e-024, 6.55197242e-038, 1.09414605e-054,
                                       1.40343441e-074, 1.38268273e-097, 1.04632487e-123, 6.08168818e-153])

    assert np.allclose(y1_double_fit.value[::10], y1_double_fit_expected, atol=1e-5)


def tie_center(model):
    """ Dummy method for testing passing of tied parameter """
    mean = 50 * model.stddev
    return mean


def test_fixed_parameters():
    """
    Test to confirm fixed parameters do not change.
    """

    x = np.linspace(0., 10., 200)
    y = 3 * np.exp(-0.5 * (x - 6.3)**2 / 0.8**2)
    y += np.random.normal(0., 0.2, x.shape)
    spectrum = Spectrum1D(flux=y*u.Jy, spectral_axis=x*u.um)

    # Test passing fixed and bounds parameters
    g_init = models.Gaussian1D(amplitude=3.*u.Jy, mean=6.1*u.um, stddev=1.*u.um,
                               fixed={'mean': True},
                               bounds={'amplitude': (2, 5)*u.Jy})

    g_fit = fit_lines(spectrum, g_init)

    assert g_fit.mean == 6.1*u.um
    assert g_fit.bounds == g_init.bounds

    # Test passing of tied parameter
    g_init = models.Gaussian1D(amplitude=3.*u.Jy, mean=6.1*u.um, stddev=1.*u.um,
                               tied={'mean': tie_center})
    g_fit = fit_lines(spectrum, g_init)

    assert g_fit.tied == g_init.tied


def test_ignore_units():
    """
    Ignore the units
    """

    #
    #  Ignore the units based on flag
    #

    # Create the spectrum
    x_single, y_single = single_peak()
    s_single = Spectrum1D(flux=y_single*u.Jy, spectral_axis=x_single*u.um)

    # Fit the spectrum
    g_init = models.Gaussian1D(amplitude=3.*u.Jy, mean=6.1*u.um, stddev=1.*u.um)
    g_fit = fit_lines(s_single, g_init, ignore_units=True)
    y_single_fit = g_fit(x_single)

    # Comparing every 10th value.
    y_single_fit_expected = np.array([3.69669474e-13, 3.57992454e-11, 2.36719426e-09, 1.06879318e-07,
               3.29498310e-06, 6.93605383e-05, 9.96945607e-04, 9.78431032e-03,
               6.55675141e-02, 3.00017760e-01, 9.37356842e-01, 1.99969007e+00,
               2.91286375e+00, 2.89719280e+00, 1.96758892e+00, 9.12412206e-01,
               2.88900005e-01, 6.24602556e-02, 9.22061121e-03, 9.29427266e-04])

    assert np.allclose(y_single_fit[::10], y_single_fit_expected, atol=1e-5)


    #
    #  Ignore the units based on there not being units on the model
    #

    # Create the spectrum
    x_single, y_single = single_peak()
    s_single = Spectrum1D(flux=y_single*u.Jy, spectral_axis=x_single*u.um)

    # Fit the spectrum
    g_init = models.Gaussian1D(amplitude=3, mean=6.1, stddev=1.)
    g_fit = fit_lines(s_single, g_init, ignore_units=True)
    y_single_fit = g_fit(x_single)

    # Comparing every 10th value.
    y_single_fit_expected = np.array([3.69669474e-13, 3.57992454e-11, 2.36719426e-09, 1.06879318e-07,
                                      3.29498310e-06, 6.93605383e-05, 9.96945607e-04, 9.78431032e-03,
                                      6.55675141e-02, 3.00017760e-01, 9.37356842e-01, 1.99969007e+00,
                                      2.91286375e+00, 2.89719280e+00, 1.96758892e+00, 9.12412206e-01,
                                      2.88900005e-01, 6.24602556e-02, 9.22061121e-03, 9.29427266e-04])

    assert np.allclose(y_single_fit[::10], y_single_fit_expected, atol=1e-5)

    #
    # Ignore the units based on not being in the model
    #

    # Create the spectrum to fit
    x_double, y_double = double_peak()
    s_double = Spectrum1D(flux=y_double*u.Jy, spectral_axis=x_double*u.um)

    # Fit the spectrum
    g1_init = models.Gaussian1D(amplitude=2.3, mean=5.6, stddev=0.1)
    g2_init = models.Gaussian1D(amplitude=1., mean=4.4, stddev=0.1)
    g12_fit = fit_lines(s_double, g1_init+g2_init, ignore_units=True)
    y12_double_fit = g12_fit(x_double)

    # Comparing every 10th value.
    y12_double_fit_expected = np.array([2.86790780e-130, 2.12984643e-103, 1.20060032e-079, 5.13707226e-059,
                                        1.66839912e-041, 4.11292970e-027, 7.69608184e-016, 1.09308800e-007,
                                        1.17844042e-002, 9.64333366e-001, 6.04322205e-002, 2.22653307e+000,
                                        5.51964567e-005, 8.13581859e-018, 6.37320251e-038, 8.85834856e-055,
                                        1.05230522e-074, 9.48850399e-098, 6.49412764e-124, 3.37373489e-153])

    assert np.allclose(y12_double_fit[::10], y12_double_fit_expected, atol=1e-5)
