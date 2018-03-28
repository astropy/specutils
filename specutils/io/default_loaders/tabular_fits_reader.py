import logging
import os

import six
from astropy.io import fits
from astropy.nddata import StdDevUncertainty
from astropy.table import Table
from astropy.units import Unit
from astropy.wcs import WCS

from ...spectra import Spectrum1D
from ..registers import data_loader


def identify_tabular_fits(origin, *args, **kwargs):
    # check if file can be opened with this reader
    # args[0] = filename
    # fits.open(args[0]) = hdulist
    return (isinstance(args[0], six.string_types) and
            # check if file is .fits
            os.path.splitext(args[0].lower())[1] == '.fits' and
            # check hdulist has more than one extension
            len(fits.open(args[0])) > 1 and
            # check if fits has BinTable extension
            isinstance(fits.open(args[0])[1], fits.BinTableHDU)
            )

@data_loader("tabular-fits", identifier=identify_tabular_fits,
             dtype=Spectrum1D)
def tabular_fits(file_name, **kwargs):
    logging.info("Spectrum file looks like tabular-fits")
    # name is not used; what was it for?
    # name = os.path.basename(file_name.rstrip(os.sep)).rsplit('.', 1)[0]

    with fits.open(file_name, **kwargs) as hdulist:
        header = hdulist[0].header

        tab = Table.read(hdulist)

        meta = {'header': header}
        wcs = WCS(header)
        uncertainty = StdDevUncertainty(tab["err"])
        data = tab["flux"] * Unit("Jy")

    return Spectrum1D(flux=data, wcs=wcs, uncertainty=uncertainty, meta=meta)
