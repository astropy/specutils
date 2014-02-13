from .. import write_fits

import numpy as np
from astropy import units as u

def test_generate_1d_header_fromdisparray():
    arr = np.linspace(5,10)*u.um
    hd1 = write_fits.generate_1d_header_fromdisparray(arr)
    hd2 = write_fits.generate_1d_header_fromdisparray(arr, reference=6*u.um)
