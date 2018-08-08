from astropy.modeling import models, fitting
import astropy.units as u
import numpy as np

from specutils.tests.spectral_examples import simulated_spectra
from specutils.spectra import Spectrum1D, SpectralRegion
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


# def test_single_peak_fit():
#     """
#     Single Peak fit.
#     """
# 
#     # Create the spectrum
#     x_single, y_single = single_peak()
#     s_single = Spectrum1D(flux=y_single*u.Jy, spectral_axis=x_single*u.um)
# 
#     # Fit the spectrum
#     g_init = models.Gaussian1D(amplitude=3.*u.Jy, mean=6.1*u.um, stddev=1.*u.um)
#     g_fit = fit_lines(s_single, g_init)
#     y_single_fit = g_fit(x_single*u.um)
# 
#     # Comparing every 10th value.
#     y_single_fit_expected = np.array([3.69669474e-13, 3.57992454e-11, 2.36719426e-09, 1.06879318e-07,
#                3.29498310e-06, 6.93605383e-05, 9.96945607e-04, 9.78431032e-03,
#                6.55675141e-02, 3.00017760e-01, 9.37356842e-01, 1.99969007e+00,
#                2.91286375e+00, 2.89719280e+00, 1.96758892e+00, 9.12412206e-01,
#                2.88900005e-01, 6.24602556e-02, 9.22061121e-03, 9.29427266e-04]) * u.Jy
# 
#     assert np.allclose(y_single_fit.value[::10], y_single_fit_expected.value, atol=1e-5)
# 
# 
# def test_single_peak_fit_window():
#     """
#     Single Peak fit with a window specified
#     """
# 
#     # Create the sepctrum
#     x_single, y_single = single_peak()
#     s_single = Spectrum1D(flux=y_single*u.Jy, spectral_axis=x_single*u.um)
# 
#     # Fit the spectrum
#     g_init = models.Gaussian1D(amplitude=3.*u.Jy, mean=5.5*u.um, stddev=1.*u.um)
#     g_fit = fit_lines(s_single, g_init, window=2*u.um)
#     y_single_fit = g_fit(x_single*u.um)
# 
#     # Comparing every 10th value.
#     y_single_fit_expected = np.array([3.69669474e-13, 3.57992454e-11, 2.36719426e-09, 1.06879318e-07,
#                                       3.29498310e-06, 6.93605383e-05, 9.96945607e-04, 9.78431032e-03,
#                                       6.55675141e-02, 3.00017760e-01, 9.37356842e-01, 1.99969007e+00,
#                                       2.91286375e+00, 2.89719280e+00, 1.96758892e+00, 9.12412206e-01,
#                                       2.88900005e-01, 6.24602556e-02, 9.22061121e-03, 9.29427266e-04]) * u.Jy
# 
#     assert np.allclose(y_single_fit.value[::10], y_single_fit_expected.value, atol=1e-5)
# 
# 
# def test_single_peak_fit_tuple_window():
#     """
#     Single Peak fit with a window specified as a tuple
#     """
# 
#     # Create the spectrum to fit
#     x_single, y_single = single_peak()
#     s_single = Spectrum1D(flux=y_single*u.Jy, spectral_axis=x_single*u.um)
# 
#     # Fit the spectrum
#     g_init = models.Gaussian1D(amplitude=3.*u.Jy, mean=5.5*u.um, stddev=1.*u.um)
#     g_fit = fit_lines(s_single, g_init, window=(6*u.um, 7*u.um))
#     y_single_fit = g_fit(x_single*u.um)
# 
#     # Comparing every 10th value.
#     y_single_fit_expected = np.array([2.29659653e-16, 6.65481577e-14, 1.20590250e-11, 1.36651116e-09,
#                                       9.68364806e-08, 4.29130753e-06, 1.18922795e-04, 2.06094025e-03,
#                                       2.23352349e-02, 1.51370231e-01, 6.41527460e-01, 1.70025811e+00,
#                                       2.81798895e+00, 2.92071072e+00, 1.89305238e+00, 7.67293541e-01,
#                                       1.94484629e-01, 3.08271836e-02, 3.05567558e-03, 1.89411116e-04]) * u.Jy
# 
#     assert np.allclose(y_single_fit.value[::10], y_single_fit_expected.value, atol=1e-5)
# 
# 
def test_double_peak_fit():
    """
    Double Peak fit.
    """

    # Create the spectrum to fit
    x_double, y_double = double_peak()
    s_double = Spectrum1D(flux=y_double*u.Jy, spectral_axis=x_double*u.um)

    # Fit the spectrum
    g1_init = models.Gaussian1D(amplitude=2.*u.Jy, mean=5.8*u.um, stddev=0.2*u.um)
    g2_init = models.Gaussian1D(amplitude=1.*u.Jy, mean=4.5*u.um, stddev=0.2*u.um)
    g12_fit = fit_lines(s_double, g1_init+g2_init)
    y12_double_fit = g12_fit(x_double*u.um)

    # Comparing every 10th value.
    y12_double_fit_expected = np.array([3.28799137e-135, 2.25240730e-109, 2.71956422e-086, 5.78745789e-066,
                                        2.17076347e-048, 1.43507070e-033, 1.67213072e-021, 3.43402371e-012,
                                        1.24300388e-005, 7.93008453e-002, 8.97346404e-001, 1.87433540e+000,
                                        8.53522418e-004, 5.85266585e-013, 4.12475472e-024, 7.88803120e-037,
                                        2.65905236e-052, 1.57986979e-070, 1.65444309e-091, 3.05364471e-115])

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
    y2_double_fit_expected = np.array([5.68946085e-113, 1.26300793e-091, 1.57494570e-072, 1.10318828e-055,
                                       4.34068462e-041, 9.59381569e-029, 1.19110227e-018, 8.30675264e-011,
                                       3.25415359e-005, 7.16092886e-002, 8.85167001e-001, 6.14618480e-002,
                                       2.39723471e-005, 5.25218534e-011, 6.46389173e-019, 4.46861199e-029,
                                       1.73530478e-041, 3.78532718e-056, 4.63826574e-073, 3.19251079e-092])

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
    y2_double_fit_expected = np.array([3.39382364e-300, 1.08912858e-247, 3.10962348e-200, 7.89905672e-158,
                                       1.78517760e-120, 3.58943748e-088, 6.42111030e-061, 1.02195704e-038,
                                       1.44708458e-021, 1.82303063e-009, 2.04330305e-002, 2.03755985e+000,
                                       1.80770157e-003, 1.42686165e-011, 1.00201856e-024, 6.26049165e-043,
                                       3.48000513e-066, 1.72103657e-094, 7.57249652e-128, 2.96433197e-166])

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
    yl_double_fit_expected = np.array([3.16533990e-300, 1.02820429e-247, 2.96806827e-200, 7.61384057e-158,
                                       1.73567772e-120, 3.51616961e-088, 6.33003280e-061, 1.01269233e-038,
                                       1.43974143e-021, 1.81897570e-009, 2.04223028e-002, 2.03759536e+000,
                                       1.80661907e-003, 1.42347766e-011, 9.96713038e-025, 6.20190400e-043,
                                       3.42937799e-066, 1.68515980e-094, 7.35871940e-128, 2.85560946e-166])

    assert np.allclose(yl_double_fit.value[::10], yl_double_fit_expected, atol=1e-5)

    # Comparing every 10th value.
    yr_double_fit_expected = np.array([3.37828686e-300, 1.08496911e-247, 3.09989636e-200, 7.87926886e-158,
                                       1.78169603e-120, 3.58418365e-088, 6.41439710e-061, 1.02124602e-038,
                                       1.44648476e-021, 1.82266264e-009, 2.04318490e-002, 2.03759561e+000,
                                       1.80774542e-003, 1.42680785e-011, 1.00184990e-024, 6.25819044e-043,
                                       3.47779387e-066, 1.71936407e-094, 7.56207196e-128, 2.95884838e-166])

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
    y1_double_fit_expected = np.array([0.00000000e+000, 0.00000000e+000, 4.60034434e-275, 6.45523250e-217,
                                       1.12263985e-165, 2.41978088e-121, 6.46425744e-084, 2.14026853e-053,
                                       8.78263645e-030, 4.46671830e-013, 2.81552405e-003, 2.19956568e+000,
                                       2.12971578e-004, 2.55572064e-015, 3.80112481e-033, 7.00677057e-058,
                                       1.60077644e-089, 4.53263229e-128, 1.59065982e-173, 6.91848736e-226])

    assert np.allclose(y1_double_fit.value[::10], y1_double_fit_expected, atol=1e-5)

    # Comparing every 10th value.
    y2_double_fit_expected = np.array([1.67318992e-125, 9.01889207e-102, 1.53067297e-080, 8.17961590e-062,
                                       1.37627329e-045, 7.29117838e-032, 1.21621972e-020, 6.38774585e-012,
                                       1.05634186e-005, 5.50025248e-002, 9.01741451e-001, 4.65481810e-002,
                                       7.56562252e-006, 3.87175525e-012, 6.23867549e-021, 3.16517823e-032,
                                       5.05621131e-046, 2.54315913e-062, 4.02757028e-081, 2.00832499e-102])

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
    y1_double_fit_expected = np.array([1.97762879e-91, 3.20701123e-74, 8.04528668e-59, 3.12225166e-45,
                                       1.87447487e-33, 1.74091144e-23, 2.50126135e-15, 5.55938567e-09,
                                       1.91152452e-04, 1.01675959e-01, 8.36646547e-01, 1.06500367e-01,
                                       2.09722740e-04, 6.38888809e-09, 3.01085871e-15, 2.19503176e-23,
                                       2.47557779e-33, 4.31914402e-45, 1.16574653e-58, 4.86738741e-74])


    assert np.allclose(y1_double_fit.value[::10], y1_double_fit_expected, atol=1e-5)
