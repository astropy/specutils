import os
import numpy as np
from specutils.io.read_IRAF_spec import make_IRAF_wave
import astropy.io.fits as fits
from specutils.io import read_fits

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_multispec_chebyshev():
    filename = data_path("spline.fits")
    spectra = read_fits.read_fits_multispec_to_list(filename)
    spec = spectra[0]
    hdus = fits.open(filename)
    flux = hdus[0].data
    meta = hdus[0].header
    wave = make_IRAF_wave(meta)[0]
    np.testing.assert_allclose(wave, spec.dispersion.value)

def test_linear_spline():
    filename = data_path("spline.fits")
    hdus = fits.open(filename)
    flux = hdus[0].data
    meta = hdus[0].header
    wave = make_IRAF_wave(meta)[0]
    c = np.array([8.61795084,  5.18731657,  0.78303947,  8.60234236,  9.65883078,
         9.86579257,  8.2952905,  1.47467708,  7.86243508,  4.54125793,
         7.96101273])
    npieces = 10
    pmin = 1
    pmax = 2304
    s = np.linspace(0, npieces, pmax-pmin+1)
    j = np.array(map(int, s))
    a = (j+1) - s
    b = s - j
    w = c[j] * a + np.take(c, 1 + j, mode='clip') * b
    hdus.close()
    np.testing.assert_allclose(wave, w)