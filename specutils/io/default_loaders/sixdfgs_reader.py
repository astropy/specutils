import astropy.io.fits as fits
from astropy.table import Table
from astropy.units import Quantity, Unit
from specutils.io.registers import data_loader
from specutils import Spectrum1D

SIXDFGS_PRIMARY_HDU_KEYWORDS = ["OBSRA", "OBSDEC", "Z", "Z_HELIO", "QUALITY"]
COUNTS_PER_SECOND = Unit("counts/s", parse_strict="silent")


def identify_6dfgs_tabular_fits(origin, *args, **kwargs):
    """
    Identify if the current file is a 6dFGS file (stored as a table)
    """
    with fits.open(args[0]) as hdulist:
        if len(hdulist) != 2:
            return False
        primary_hdu = hdulist[0]
        if primary_hdu.header["NAXIS"] != 0:
            return False
        for keyword in SIXDFGS_PRIMARY_HDU_KEYWORDS:
            if keyword not in primary_hdu.header:
                return False
        return True


@data_loader("6dFGS-tabular",
             identifier=identify_6dfgs_tabular_fits, dtype=Spectrum1D,
             extensions=["fit", "fits"])
def sixdfgs_tabular_fits_loader(file_obj, **kwargs):
    """
    Load the tabular variant of a 6dF Galaxy Survey (6dFGS) file.

    6dFGS used the Six-degree Field instrument on the UK Schmidt Telescope
    (UKST) at the Siding Spring Observatory (SSO) near Coonabarabran,
    Australia. Further details can be found at http://www.6dfgs.net/, or
    https://docs.datacentral.org.au/6dfgs/. Catalogues and spectra were
    produced, with the spectra being provided as both fits tables and as fits
    images. This loads the tabular variant of the spectra. Note that this does
    not include uncertainties - uncertainties are only provided in the image
    variants.

    Parameters
    ----------
    file_obj: str, file-like or HDUList
         FITS file name, object (provided from name by Astropy I/O Registry),
         or HDUList (as resulting from astropy.io.fits.open()).

    Returns
    -------
    data: Spectrum1D
        The 6dF spectrum that is represented by the data in this table.
    """
    if isinstance(file_obj, fits.hdu.hdulist.HDUList):
        hdulist = file_obj
    else:
        hdulist = fits.open(file_obj, **kwargs)

    header = hdulist[0].header
    table = Table.read(hdulist)
    flux = Quantity(table["FLUX"])
    wavelength = Quantity(table["WAVE"])

    if flux.unit == COUNTS_PER_SECOND:
        flux._unit = Unit("count/s")
    meta = {"header": header}

    if not isinstance(file_obj, fits.hdu.hdulist.HDUList):
        hdulist.close()

    return Spectrum1D(flux=flux, spectral_axis=wavelength, meta=meta)
