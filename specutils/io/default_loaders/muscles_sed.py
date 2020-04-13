import logging
import os

from astropy.io import fits
from astropy.nddata import StdDevUncertainty
from astropy.table import Table
from astropy.units import Unit, Quantity
from astropy.wcs import WCS

from ...spectra import Spectrum1D
from ..registers import data_loader


def identify_muscles_sed(origin, *args, **kwargs):
    # check if file can be opened with this reader
    # args[0] = filename
    # fits.open(args[0]) = hdulist
    with fits.open(args[0]) as hdulist:
        # Test if fits has extension of type BinTable and check against
        # known keys of already defined specific formats
        return (len(hdulist) > 1 and
                isinstance(hdulist[1], fits.BinTableHDU) and
                fits.getheader(args[0]).get('TELESCOP') == 'MULTI' and
                fits.getheader(args[0]).get('HLSPACRN') == 'MUSCLES' and
                fits.getheader(args[0]).get('PROPOSID') == 13650)


@data_loader(label="MUSCLES SED", identifier=identify_muscles_sed,
             dtype=Spectrum1D, extensions=['fits'])
def muscles_sed(file_obj, **kwargs):
    """
    Load spectrum from a MUSCLES Treasury Survey panchromatic SED FITS file.

    Parameters
    ----------
    file_obj: str, file-like, or HDUList
          FITS file name, object (provided from name by Astropy I/O Registry),
          or HDUList (as resulting from astropy.io.fits.open()).

    Returns
    -------
    data: Spectrum1D
        The spectrum that is represented by the data in this table.
    """
    # name is not used; what was it for?
    # name = os.path.basename(file_name.rstrip(os.sep)).rsplit('.', 1)[0]
    if isinstance(file_obj, fits.hdu.hdulist.HDUList):
        hdulist = file_obj
    else:
        hdulist = fits.open(file_obj, **kwargs)

    header = hdulist[0].header

    tab = Table.read(hdulist[1])

    meta = {'header': header}
    uncertainty = StdDevUncertainty(tab["ERROR"])
    data = Quantity(tab["FLUX"])
    wavelength = Quantity(tab["WAVELENGTH"])

    if not isinstance(file_obj, fits.hdu.hdulist.HDUList):
        hdulist.close()

    return Spectrum1D(flux=data, spectral_axis=wavelength,
                      uncertainty=uncertainty, meta=meta)
