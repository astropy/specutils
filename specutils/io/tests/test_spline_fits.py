import os
import numpy as np
from specutils.io import read_fits
import pytest


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_multispec_linear_spline():
    pytest.importorskip("scipy")
    filename = data_path("spline.fits")
    spectra = read_fits.read_fits_spectrum1d(filename)
    c1 = np.array([8.61795084, 5.18731657, 0.78303947, 8.60234236,
                   9.65883078, 9.86579257, 8.2952905, 1.47467708,
                   7.86243508, 4.54125793, 7.96101273])
    npieces = 10
    pmin = 1
    pmax = 2304
    pixels = np.arange(pmin, pmax + 1) * 1.0
    s = (pixels - pmin) / (pmax - pmin) * npieces
    j = s.astype(int)
    a = (j + 1) - s
    b = s - j
    w1 = np.take(c1, j, mode='clip') * a + np.take(c1, 1 + j, mode='clip') * b
    np.testing.assert_allclose(w1, spectra[0].dispersion.value)


# expected to fail as the formula given in http://iraf.net/irafdocs/specwcs.php
# is not matching the result from scipy's cubic spline implementation
@pytest.mark.xfail
def test_multispec_cubic_spline():
    pytest.importorskip("scipy")
    filename = data_path("spline.fits")
    spectra = read_fits.read_fits_spectrum1d(filename)
    npieces = 10
    pmin = 1
    pmax = 2304
    pixels = np.arange(pmin, pmax + 1) * 1.0
    s = (pixels - pmin) / (pmax - pmin) * npieces
    j = np.array(map(int, s))
    a = (j + 1) - s
    b = s - j

    c2 = np.array([1.91932706, 4.50545551, 8.47358561, 8.37680956,
                   3.0511021, 2.54651656, 7.5758634, 7.68131867, 7.58694718,
                   3.14098023, 7.70766882, 7.90089733, 9.80179082])

    x0 = a ** 3
    x1 = (1 + 3 * a * (1 + a * b))
    x2 = (1 + 3 * b * (1 + a * b))
    x3 = b ** 3
    w2 = np.take(c2, j, mode='clip') * x0 + np.take(c2, 1 + j, mode='clip') * x1 \
        + np.take(c2, 2 + j, mode='clip') * x2 + np.take(c2, 3 + j,
                                                         mode='clip') * x3
    np.testing.assert_allclose(w2, spectra[1].dispersion.value)
