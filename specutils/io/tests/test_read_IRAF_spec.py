import astropy.io.fits as fits
import astropy.io.ascii as ascii
from ..read_IRAF_spec import read_IRAF_spec

from astropy.io.ascii.tests.common import  assert_almost_equal



def test_multispec_legendre():
    iraf = ascii.read('t/TRES.dat', data_start = 127, Reader = ascii.NoHeader, names = ['wave', 'flux'])
    spec = read_IRAF_spec('t/TRES.fits')
    assert_almost_equal(iraf['wave'], spec.dispersion[10,:])

def test_multispec_chebyshev():
    iraf = ascii.read('t/AAO_11.txt', data_start = 175, Reader = ascii.NoHeader, names = ['wave', 'flux'])
    spec = read_IRAF_spec('t/AAO.fits')
    assert_almost_equal(iraf['wave'], spec.dispersion[10,:])
