from astropy.io import fits
from astropy.nddata import StdDevUncertainty
from astropy.table import Table
from astropy.units import Quantity

from ...spectra import Spectrum1D
from ..registers import data_loader
from ..parsing_utils import read_fileobj_or_hdulist


def identify_muscles_sed(origin, *args, **kwargs):
    # check if file can be opened with this reader
    # args[0] = filename
    # fits.open(args[0]) = hdulist
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        # Test if fits has extension of type BinTable and check against
        # known keys of already defined specific formats
        return (len(hdulist) > 1 and
                isinstance(hdulist[1], fits.BinTableHDU) and
                hdulist[0].header.get('TELESCOP') == 'MULTI' and
                hdulist[0].header.get('HLSPACRN') == 'MUSCLES' and
                hdulist[0].header.get('PROPOSID') == 13650)


@data_loader(
    label="MUSCLES SED", identifier=identify_muscles_sed, dtype=Spectrum1D,
    extensions=['fits'], priority=10,
)
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

    with read_fileobj_or_hdulist(file_obj, **kwargs) as hdulist:
        header = hdulist[0].header

        tab = Table.read(hdulist[1])

    meta = {'header': header}
    uncertainty = StdDevUncertainty(tab["ERROR"])
    data = Quantity(tab["FLUX"])
    wavelength = Quantity(tab["WAVELENGTH"])

    return Spectrum1D(flux=data, spectral_axis=wavelength,
                      uncertainty=uncertainty, meta=meta)
