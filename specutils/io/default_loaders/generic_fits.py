import os

import six
from astropy.io import fits
from astropy.nddata import StdDevUncertainty
from astropy.table import Table
from astropy.units import Unit
from astropy.wcs import WCS

from specutils.io.registers import data_loader
from specutils.spectra import Spectrum1D


def identify_generic_fits(origin, *args, **kwargs):
    return (isinstance(args[0], six.string_types) and
            os.path.splitext(args[0].lower())[1] == '.fits')


@data_loader("generic fits", identifier=identify_generic_fits,
             dtype=Spectrum1D)
def generic_fits(file_name, **kwargs):
    name = os.path.basename(file_name.rstrip(os.sep)).rsplit('.', 1)[0]

    with fits.open(file_name, **kwargs) as hdulist:
        header = hdulist[0].header

        tab = Table.read(file_name)

        meta = {'header': header}
        wcs = WCS(hdulist[0].header)
        uncertainty = StdDevUncertainty(tab["err"])
        data = tab["flux"] * Unit("Jy")

    return Spectrum1D(flux=data, wcs=wcs, uncertainty=uncertainty, meta=meta)