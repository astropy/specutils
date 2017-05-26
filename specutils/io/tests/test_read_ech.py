import astropy.io.fits as fits
import astropy.io.ascii as ascii
from rdech import rdech

from astropy.io.ascii.tests.common import  assert_almost_equal

def test_read_ech():
    '''Compare with a IDL REDUCE output.
    
    This output is saved here as a txt file, so that you don't need the IDL code
    to run the comparison
    
    Below is the idl code used to make the comparison file (assuming that
    REDUCE is in your idl path)::
    
        idl> ; input
        idl> file = 't/REDUCE.fits.ech'
        idl> rdech, ech, file
        idl>
        idl> ; reformat for output
        idl> out = make_array(3,2048)
        idl> out[0,*] = ech.WAVE[*,10]
        idl> out[1,*] = ech.SPEC[*,10]
        idl> out[2,*] = ech.SIG[*,10]
        idl> print, out , FORMAT='(F7.2,1X,F7.5,1X,e10.5)'
        idl> 
        idl> ; write
        idl> fname='t/REDUCE.dat'
        idl> OPENW,1,fname
        idl> PRINTF, 1, out , FORMAT='(F7.2,1X,F7.5,1X,e10.5)'
        idl> CLOSE,1
        
    This file has no RADVEL correction, there is only a header entry for
    the barycentric correction.
    '''
    reduce_orig = ascii.read('t/REDUCE.dat', Reader = ascii.NoHeader, names = ['wave', 'flux', 'sigma'])
    spec = rdech('t/REDUCE.fits.ech')
    
    assert_almost_equal(spec.dispersion[:,10], reduce_orig['wave'])
    assert_almost_equal(spec.data[:,10], reduce_orig['flux'])
    assert_almost_equal(spec.uncertainty.array[:,10], reduce_orig['sigma'])
