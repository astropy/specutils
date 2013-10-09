import os
import astropy.io.ascii as ascii
import numpy as np
from ..read_IRAF_spec import read_IRAF_spec

from astropy.io.ascii.tests.common import  assert_almost_equal


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 't')
    return os.path.join(data_dir, filename)

def test_multispec_legendre():
    iraf = ascii.read(data_path('TRES.dat'), data_start = 127, Reader = ascii.NoHeader, names = ['wave', 'flux'])
    spec = read_IRAF_spec(data_path('TRES.fits'))
    assert_almost_equal(iraf['wave'], spec.dispersion[10,:])

def test_multispec_chebyshev():
    iraf = ascii.read(data_path('AAO_11.txt'), data_start = 175, Reader = ascii.NoHeader, names = ['wave', 'flux'])
    spec = read_IRAF_spec(data_path('AAO.fits'))
    assert_almost_equal(iraf['wave'], spec.dispersion[10,:])

def test_1dspec_UVES():
    spec = read_IRAF_spec(data_path('UVES.fits'))
    assert np.abs(spec.dispersion[0,:] - 3732.056) < 0.01

def test_1dspec_vrad():
    spec = read_IRAF_spec(data_path('1d.fits'))
    assert np.abs(spec.dispersion[0,:] - 509.46) < 0.1
    assert spec.unit == 'km/s'
