from astropy.units import Unit
from astropy.nddata import StdDevUncertainty

from ...spectra import Spectrum
from ..registers import data_loader
from ..parsing_utils import read_fileobj_or_hdulist

__all__ = ['cos_identify', 'cos_spectrum_loader']


def cos_identify(origin, *args, **kwargs):
    """Check whether given file contains HST/COS spectral data."""
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        return (hdulist[0].header['TELESCOP'] == 'HST' and hdulist[0].header['INSTRUME'] == 'COS')


@data_loader(
    label="HST/COS", identifier=cos_identify,
    extensions=['FITS', 'FIT', 'fits', 'fit'], priority=10,
)
def cos_spectrum_loader(file_obj, **kwargs):
    """
    Load COS spectral data from the MAST archive into a spectrum object.

    Parameters
    ----------
    file_obj: str, file-like, or HDUList
         FITS file name, object (provided from name by Astropy I/O Registry),
         or HDUList (as resulting from astropy.io.fits.open()).

    Returns
    -------
    data: Spectrum
        The spectrum that is represented by the data in this table.
    """

    with read_fileobj_or_hdulist(file_obj, **kwargs) as hdulist:
        header = hdulist[0].header
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

    return Spectrum(flux=data,
                      spectral_axis=dispersion,
                      uncertainty=uncertainty,
                      meta=meta)
