from ..IUE import read_IUE_mxlo
from ..IUE import ApertureException

#already written for astropy, if might not quite fit here
from astropy.io.ascii.tests.common import  assert_almost_equal
import pytest

def test_mxlo():
    spec1 = read_IUE_mxlo('t/lwp11854.mxlo')
    assert spec1.meta['LDATEOBS'] == '11/10/87'

    
    spec3 = read_IUE_mxlo('t/swp02283.mxlo', aperture = 'SMALL')
    spec4 = read_IUE_mxlo('t/swp02283.mxlo', aperture = 'LARGE')
    assert spec3.dispersion[0] == 1050.
    assert abs(spec3.dispersion[1]- (1050.+1.676338)) < 1e-4


def test_multiple_apertures_exception():
    pytest.raises(ApertureException, read_IUE_mxlo, 't/swp02283.mxlo')
    pytest.raises(ApertureException, read_IUE_mxlo, 't/swp02283.mxlo', aperture = 'XXX')
