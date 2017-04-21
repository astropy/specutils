import os

import six
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS

from specutils.io.registers import data_loader
from specutils.spectra import Spectrum1D


def identify_wcs1d_fits(origin, *args, **kwargs):
    return (isinstance(args[0], six.string_types) and
            os.path.splitext(args[0].lower())[1] == '.fits' and
            fits.getheader(args[0])['NAXIS'] == 1 and
            'CTYPE1' in fits.getheader(args[0])
           )


@data_loader("wcs1d-fits-reader", identifier=identify_wcs1d_fits,
             dtype=Spectrum1D)
def wcs1d_fits(file_name, **kwargs):
    # name is not used; what was it for?
    # name = os.path.basename(file_name.rstrip(os.sep)).rsplit('.', 1)[0]

    with fits.open(file_name, **kwargs) as hdulist:
        header = hdulist[0].header
        wcs = WCS(header)

        if 'BUNIT' in header:
            data = u.Quantity(hdulist[0].data, unit=header['BUNIT'])
        else:
            data = hdulist[0].data

        meta = {'header': header}

    return Spectrum1D(flux=data, wcs=wcs, meta=meta)

