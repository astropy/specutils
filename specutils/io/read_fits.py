# a simple fits readers that puts a lot of emphasis on the specwcs

from astropy.io import fits

from specutils.wcs import specwcs
from specutils import Spectrum1D

def read_fits(filename):
    data = fits.open(filename)
    header = fits.getheader(filename)

    for fits_wcs in specwcs.fits_capable_wcs:
        try:
            wcs = fits_wcs.from_fits_header(header)
        except specwcs.Spectrum1DWCSError:
            continue
        else:
            break
    else:
        raise specwcs.Spectrum1DWCSError('File %s does not contain a spectutils-readable WCS - exiting')

    return Spectrum1D(data, wcs=wcs)



