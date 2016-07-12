from specviz.interfaces.decorators import data_loader
from specviz.core.data import Spectrum1DRef

import os

from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.units import Unit
from astropy.nddata import StdDevUncertainty


def fits_identify(*args, **kwargs):
    """
    Check whether given filename is FITS. This is used for Astropy I/O
    Registry.
    """
    return (isinstance(args[0], str) and
            args[0].lower().split('.')[-1] in ['fits', 'fit'])


@data_loader(label="Simple Fits", identifier=fits_identify)
def simple_generic_loader(file_name, **kwargs):
    name = os.path.basename(file_name.name.rstrip(os.sep)).rsplit('.', 1)[0]
    hdulist = fits.open(file_name, **kwargs)

    header = hdulist[0].header

    tab = Table.read(file_name)

    meta = {'header': header}
    wcs = WCS(hdulist[0].header)
    unit = Unit('Jy')
    uncertainty = StdDevUncertainty(tab["err"])
    data = tab["flux"]

    hdulist.close()

    return Spectrum1DRef(data=data, name=name, wcs=wcs,
                         uncertainty=uncertainty, unit=unit, meta=meta)

