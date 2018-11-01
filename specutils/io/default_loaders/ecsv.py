import os

from astropy.table import Table

from specutils.io.registers import data_loader, custom_writer
from specutils import Spectrum1D

__all__ = ['ecsv_identify', 'ecsv_spectrum_loader']


def ecsv_identify(*args, **kwargs):
    """Check if it's an ECSV file."""
    name = os.path.basename(args[0])

    if name.lower().split('.')[-1] == 'ecsv':
       return True

    return False


@data_loader(label="Spectrum ECSV", identifier=ecsv_identify, extensions=['ecsv'])
def ecsv_spectrum_loader(file_name, **kwargs):
    """
    Load spectrum from ECSV file

    Parameters
    ----------
    file_name: str
        The path to the ECSV file

    Returns
    -------
    data: Spectrum1D
        The data.
    """
    table = Table.read(file_name, format='ascii.ecsv')

    unit = table['Intensity'].unit
    disp_unit = table['Wavelength'].unit
    flux = table['Intensity'] * unit
    dispersion = table['Wavelength'] * disp_unit

    return Spectrum1D(flux=flux,
                      spectral_axis=dispersion,
                      meta=table.meta)
