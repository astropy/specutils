import os
from specutils.io import read_fits, write_fits
import numpy as np

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_linear_write1():
    spec = read_fits.read_fits_spectrum1d(data_path('UVES.fits'))
    data = spec.data
    wcs = spec.wcs

    write_fits.write(data_path('test.fits'), data, wcs)
    test_spec = read_fits.read_fits_spectrum1d(data_path('test.fits'))
    np.testing.assert_allclose(spec.data, test_spec.data)
    np.testing.assert_allclose(spec.dispersion, test_spec.dispersion)

def test_linear_write2():
    spec = read_fits.read_fits_spectrum1d(data_path('gbt_1d.fits'))
    data = spec.data
    wcs = spec.wcs

    write_fits.write(data_path('test.fits'), data, wcs)
    test_spec = read_fits.read_fits_spectrum1d(data_path('test.fits'))
    np.testing.assert_allclose(spec.data, test_spec.data)
    np.testing.assert_allclose(spec.dispersion, test_spec.dispersion)

# TODO: create tests that utilize all methods to write fits header