import logging
import os

from astropy.io import fits
from astropy.nddata import StdDevUncertainty
from astropy.table import Table
from astropy.units import Unit
from astropy.wcs import WCS

from ...spectra import Spectrum1D
from ..registers import data_loader


def identify_muscles_sed(origin, *args, **kwargs):
    # check if file can be opened with this reader
    # args[0] = filename
    # fits.open(args[0]) = hdulist
    return (isinstance(args[0], str) and
            # check if file is .fits
            args[0].endswith('sed.fits') and
            # check hdulist has more than one extension
            len(fits.open(args[0])) > 1 and
            # check if fits has BinTable extension
            isinstance(fits.open(args[0])[1], fits.BinTableHDU) and
            # check if MUSCLES proposal ID is in fits header
            fits.open(args[0])[0].header['PROPOSID'] == 13650
            )


@data_loader("muscles-sed", identifier=identify_muscles_sed,
             dtype=Spectrum1D, extensions=['fits'])
def muscles_sed(file_name, **kwargs):
    """
    Load spectrum from a MUSCLES SED FITS file.

    Parameters
    ----------
    file_name: str
        The path to the FITS file.

    Returns
    -------
    data: Spectrum1D
        The spectrum that is represented by the data in this table.
    """
    # name is not used; what was it for?
    # name = os.path.basename(file_name.rstrip(os.sep)).rsplit('.', 1)[0]

    with fits.open(file_name, **kwargs) as hdulist:
        header = hdulist[0].header

        tab = Table.read(hdulist)

        meta = {'header': header}
        uncertainty = StdDevUncertainty(tab["ERROR"])
        data = tab["FLUX"]
        wavelength = tab["WAVELENGTH"]

    return Spectrum1D(flux=data, spectral_axis=wavelength,
                      uncertainty=uncertainty, meta=meta,
                      unit=data.unit,
                      spectral_axis_unit=wavelength.unit)
