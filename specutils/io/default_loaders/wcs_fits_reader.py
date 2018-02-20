import logging
import os

import six
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS

from ...spectra import Spectrum1D
from ..registers import data_loader


def identify_wcs1d_fits(origin, *args, **kwargs):
    # check if file can be opened with this reader
    # args[0] = filename
    return (isinstance(args[0], six.string_types) and
            os.path.splitext(args[0].lower())[1] == '.fits' and
            # check if number of axes is one
            fits.getheader(args[0])['NAXIS'] == 1 and
            # check if CTYPE1 kep is in the header
            'CTYPE1' in fits.getheader(args[0])
            )


@data_loader("wcs1d-fits", identifier=identify_wcs1d_fits,
             dtype=Spectrum1D)
def wcs1d_fits(file_name, spectral_axis_unit=None, **kwargs):
    """
       Parameters
       ----------
       file_name : str

        spectral_axis_unit: str or unit, optional
            Optional string or unit object to specify units of spectral axis.
    """
    logging.info("Spectrum file looks like wcs1d-fits")

    with fits.open(file_name, **kwargs) as hdulist:
        header = hdulist[0].header
        wcs = WCS(header)

        if 'BUNIT' in header:
            data = u.Quantity(hdulist[0].data, unit=header['BUNIT'])
        else:
            data = hdulist[0].data

        if spectral_axis_unit is not None:
            wcs.wcs.cunit[0] = spectral_axis_unit

        meta = {'header': header}

    return Spectrum1D(flux=data, wcs=wcs, meta=meta)
