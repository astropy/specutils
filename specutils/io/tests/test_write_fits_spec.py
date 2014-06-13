import os
from specutils.io import read_fits, write_fits
import numpy as np

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_linear_write():
    spec = read_fits.read_fits_spectrum1d(data_path('UVES.fits'))
    data = spec.data
    wcs = spec.wcs

    write_fits.write(data_path('test.fits'), data, wcs)
    test_spec = read_fits.read_fits_spectrum1d(data_path('test.fits'))
    np.testing.assert_allclose(spec.data, test_spec.data)
    np.testing.assert_allclose(spec.dispersion, test_spec.dispersion)

