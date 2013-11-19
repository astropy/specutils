#already written for astropy, if might not quite fit here
import pytest
pytestmark = pytest.mark.skipif(True, reason='tests need to be reimplemented')


import os
#from ..IUE import read_IUE_mxlo
#from ..IUE import ApertureException



#marking whole file for skipping (as most of this will be re-implemented using the WCS classes


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 't')
    return os.path.join(data_dir, filename)

def test_mxlo():
    spec1 = read_IUE_mxlo(data_path('lwp11854.mxlo'))
    assert spec1.meta['LDATEOBS'] == '11/10/87'

    
    spec3 = read_IUE_mxlo(data_path('swp02283.mxlo'), aperture='SMALL')
    spec4 = read_IUE_mxlo(data_path('swp02283.mxlo'), aperture='LARGE')
    assert spec3.dispersion[0].value == 1050.
    assert abs(spec3.dispersion[1]- (1050.+1.676338)) < 1e-4


def test_multiple_apertures_exception():
    pytest.raises(ApertureException, read_IUE_mxlo, data_path('swp02283.mxlo'))
    pytest.raises(ApertureException, read_IUE_mxlo, data_path('swp02283.mxlo'), aperture = 'XXX')
