import os

from astropy.io import fits
from astropy.units import Unit
from astropy.nddata import StdDevUncertainty

from specutils.io.registers import data_loader
from specutils import Spectrum1D

__all__ = ['stis_identify', 'stis_spectrum_loader']


def stis_identify(origin, *args, **kwargs):
    """Check whether given file contains HST/STIS spectral data."""
    with fits.open(args[0]) as hdulist:
        if hdulist[0].header['TELESCOP'] == 'HST' and hdulist[0].header['INSTRUME'] == 'STIS':
            return True

    return False


@data_loader(label="HST/STIS", identifier=stis_identify, extensions=['FITS', 'FIT', 'fits', 'fit'])
def stis_spectrum_loader(file_obj, **kwargs):
    """
    Load STIS spectral data from the MAST archive into a spectrum object.

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
    if isinstance(file_obj, fits.hdu.hdulist.HDUList):
        hdulist = file_obj
    else:
        hdulist = fits.open(file_obj, **kwargs)

    header = hdulist[0].header
    name = header.get('FILENAME')
    meta = {'header': header}

    unit = Unit("erg/cm**2 Angstrom s")
    disp_unit = Unit('Angstrom')
    data = hdulist[1].data['FLUX'].flatten() * unit
    dispersion = hdulist[1].data['wavelength'].flatten() * disp_unit
    uncertainty = StdDevUncertainty(hdulist[1].data["ERROR"].flatten() * unit)

    sort_idx = dispersion.argsort()
    dispersion = dispersion[sort_idx]
    data = data[sort_idx]
    uncertainty = uncertainty[sort_idx]

    if not isinstance(file_obj, fits.hdu.hdulist.HDUList):
        hdulist.close()

    return Spectrum1D(flux=data,
                      spectral_axis=dispersion,
                      uncertainty=uncertainty,
                      meta=meta)
