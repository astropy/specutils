import os

from astropy.io import fits
from astropy.units import Unit
from astropy.nddata import StdDevUncertainty

from ...interfaces import data_loader
from ...core.data import Spectrum1DRef

__all__ = ['cos_identify', 'cos_spectrum_loader']

def cos_identify(*args, **kwargs):
    """ Check whether given file contains HST/COS spectral data.

    """

    with fits.open(args[0]) as hdu:
        if hdu[0].header['TELESCOP'] == 'HST' and hdu[0].header['INSTRUME'] == 'COS':
           return True

    return False


@data_loader(label="HST/COS", priority=10, identifier=cos_identify)
def cos_spectrum_loader(file_name, **kwargs):
    """ Load file from COS spectral data into a spectrum object

    """

    name = os.path.basename(file_name)

    with fits.open(file_name, **kwargs) as hdu:
        header = hdu[0].header
        meta = {'header': header}

        uncertainty = StdDevUncertainty(hdu[1].data["ERROR"].flatten())
        data = hdu[1].data['FLUX'].flatten()
        unit = Unit("erg/cm**2 Angstrom s")

    return Spectrum1DRef(data=data, name=name,
                         uncertainty=uncertainty, unit=unit, meta=meta)
