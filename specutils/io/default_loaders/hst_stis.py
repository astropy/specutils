from astropy.units import Unit
from astropy.nddata import StdDevUncertainty

from ...spectra import Spectrum1D
from ..registers import data_loader
from ..parsing_utils import read_fileobj_or_hdulist

__all__ = ['stis_identify', 'stis_spectrum_loader']


def stis_identify(origin, *args, **kwargs):
    """Check whether given file contains HST/STIS spectral data."""
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        return (hdulist[0].header['TELESCOP'] == 'HST' and hdulist[0].header['INSTRUME'] == 'STIS')

    return False


@data_loader(
    label="HST/STIS", identifier=stis_identify,
    extensions=['FITS', 'FIT', 'fits', 'fit'], priority=10,
)
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

    return Spectrum1D(flux=data,
                      spectral_axis=dispersion,
                      uncertainty=uncertainty,
                      meta=meta)
