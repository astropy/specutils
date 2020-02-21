from astropy.io import fits
from astropy.units import Unit

from specutils.io.registers import data_loader
from specutils import XraySpectrum1D

__all__ = ['chandra_hetg_identify', 'hetg_spectrum_loader']

def chandra_hetg_identify(origin, *args, **kwargs):
    """Check whether given file contains Chandra HETG spectral data."""
    with fits.open(args[0]) as hdu:
        # Check that it is a Chandra spectrum file
        if hdu[0].header['TELESCOP'] == 'CHANDRA' and 'SPECTRUM' in [h.name for h in hdu]:
            # Check that it's an HETG file
            if hdu['SPECTRUM'].header['GRATING'] == 'HETG':
                return True
    return False

@data_loader(label="chandra_hetg", identifier=chandra_hetg_identify, extensions=['pha','pha2', 'PHA', 'PHA2', 'FITS', 'FIT', 'fits', 'fit'])
def hetg_spectrum_loader(file_name, arf=None, rmf=None):
    """
    Load Chandra HETG spectral data from a file into a spectrum object.

    Parameters
    ----------
    file_name: str
        The path to the FITS file

    arf : str OR AreaResponse
        Filename for the area response file (ARF) or a pre-loaded AreaResponse object

    rmf : str OR ResponseMatrix
        Filename for the response matrix file (RMF) or a pre-loaded ResponseMatrix object

    Returns
    -------
    data: XraySpectrum1D
        The spectrum that is represented by the data in this table.
    """

    with fits.open(file_name) as hdu:
        header = hdu[0].header
        meta   = {'header': header}
        data   = hdu[1].data

        bin_unit = Unit(data.columns['BIN_LO'].unit)
        bin_lo   = data['BIN_LO'] * bin_unit
        bin_hi   = data['BIN_HI'] * bin_unit

        counts   = data['COUNTS'] * Unit('ct')
        exposure = hdu[1].header['EXPOSURE'] * Unit('second')

    return XraySpectrum1D(bin_lo, bin_hi, counts,
                          exposure=exposure, arf=arf, rmf=rmf)
