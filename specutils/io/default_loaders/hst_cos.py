import os

from astropy.io import fits
from astropy.units import Unit
from astropy.nddata import StdDevUncertainty

from specutils.io.registers import data_loader
from specutils import Spectrum1D

__all__ = ['cos_identify', 'cos_spectrum_loader']


def cos_identify(origin, *args, **kwargs):
    """Check whether given file contains HST/COS spectral data."""
    with fits.open(args[0]) as hdu:
        if hdu[0].header['TELESCOP'] == 'HST' and hdu[0].header['INSTRUME'] == 'COS':
            return True

    return False


@data_loader(label="HST/COS", identifier=cos_identify, extensions=['FITS', 'FIT', 'fits', 'fit'])
def cos_spectrum_loader(file_name, **kwargs):
    """
    Load COS spectral data from the MAST archive into a spectrum object.

    Parameters
    ----------
    file_name: str
        The path to the FITS file

    Returns
    -------
    data: Spectrum1D
        The spectrum that is represented by the data in this table.
    """

    name = os.path.basename(file_name)

    with fits.open(file_name, **kwargs) as hdu:
        header = hdu[0].header
        meta = {'header': header}

        unit = Unit("erg/cm**2 Angstrom s")
        disp_unit = Unit('Angstrom')
        data = hdu[1].data['FLUX'].flatten() * unit
        dispersion = hdu[1].data['wavelength'].flatten() * disp_unit
        uncertainty = StdDevUncertainty(hdu[1].data["ERROR"].flatten() * unit)

        sort_idx = dispersion.argsort()
        dispersion = dispersion[sort_idx]
        data = data[sort_idx]
        uncertainty = uncertainty[sort_idx]

    return Spectrum1D(flux=data,
                      spectral_axis=dispersion,
                      uncertainty=uncertainty,
                      meta=meta)
