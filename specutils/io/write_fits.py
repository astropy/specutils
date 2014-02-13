"""
Tools for writing FITS files
As of first commit, there is no actual write mechanism, just header generators
"""
import numpy as np

from astropy.io import fits
from astropy import units as u

type_to_ctype = {'length':'WAVE',
                 'frequency':'FREQ',
                 'speed':'VELO',
                }

def generate_1d_header_fromdisparray(arr, cdelt_tolerance=1e-8, reference=None,
                                     unit=None):
    """
    Parameters
    ----------
    cdelt_tolerance : float
        Tolerance in the difference between pixels that determines
        how near to linear the dispersion axis must be
    """
    header = fits.Header()

    # convert array to numpy array with no units
    if hasattr(arr,'unit'):
        if unit is None:
            unit = arr.unit
        arr = arr.value

    # convert string to unit
    if not hasattr(unit,'to_string'):
        unit = u.Unit(unit)
    
    # determine CDELT
    dxarr = np.diff(arr)
    if abs(dxarr.max()-dxarr.min())/abs(dxarr.min()) < cdelt_tolerance:
        cdelt = dxarr.mean().flat[0]
    else:
        raise ValueError("Dispersion array is not linear.")

    # determine CRVAL, CRPIX
    crval = arr.min()
    crpix = 1

    if reference is not None:
        restfrq = reference.to(u.Hz, u.spectral())
        header['RESTFRQ'] = restfrq

    header['CRVAL1'] = crval
    header['CDELT1'] = cdelt
    header['CRPIX1'] = crpix
    header['CUNIT1'] = unit.to_string()
    header['CTYPE1'] = type_to_ctype[unit.physical_type]

    return header
