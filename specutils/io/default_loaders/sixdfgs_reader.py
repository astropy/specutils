from astropy.nddata import VarianceUncertainty
from astropy.table import Table
from astropy.units import Quantity, Unit
from astropy.wcs import WCS

from ...spectra import Spectrum1D, SpectrumList
from ..registers import data_loader
from ..parsing_utils import read_fileobj_or_hdulist

SIXDFGS_PRIMARY_HDU_KEYWORDS = ["OBSRA", "OBSDEC", "Z", "Z_HELIO", "QUALITY"]
COUNTS_PER_SECOND = Unit("counts/s", parse_strict="silent")
ANGSTROMS = Unit("Angstroms", parse_strict="silent")


def identify_6dfgs_tabular_fits(origin, *args, **kwargs):
    """
    Identify if the current file is a 6dFGS file (stored as a table)
    """
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        if len(hdulist) != 2:
            return False
        primary_hdu = hdulist[0]
        if primary_hdu.header["NAXIS"] != 0:
            return False
        for keyword in SIXDFGS_PRIMARY_HDU_KEYWORDS:
            if keyword not in primary_hdu.header:
                return False
        return True


def identify_6dfgs_split_fits(origin, *args, **kwargs):
    """
    Identify if the current file is a 6dFGS file (stored in the split variant).
    """
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        if len(hdulist) != 1:
            return False
        primary_hdu = hdulist[0]
        if primary_hdu.header["NAXIS2"] not in (3, 4):
            return False
        for keyword in SIXDFGS_PRIMARY_HDU_KEYWORDS:
            if keyword not in primary_hdu.header:
                return False
        return True


def identify_6dfgs_combined_fits(origin, *args, **kwargs):
    """
    Identify if the current file is a 6dFGS file (stored in the combined
    variant).
    """
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        if len(hdulist) < 8:
            return False
        first_spectrum_hdu = hdulist[5]
        if first_spectrum_hdu.header["NAXIS2"] not in (3, 4):
            return False
        for keyword in SIXDFGS_PRIMARY_HDU_KEYWORDS:
            if keyword not in first_spectrum_hdu.header:
                return False
        return True


@data_loader(
    "6dFGS-tabular", identifier=identify_6dfgs_tabular_fits, dtype=Spectrum1D,
    extensions=["fit", "fits"], priority=10,
)
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

    with read_fileobj_or_hdulist(file_obj, **kwargs) as hdulist:
        header = hdulist[0].header
        table = Table.read(hdulist)

    flux = Quantity(table["FLUX"])
    wavelength = Quantity(table["WAVE"])

    if flux.unit == COUNTS_PER_SECOND:
        flux._unit = Unit("count/s")
    meta = {"header": header}

    return Spectrum1D(flux=flux, spectral_axis=wavelength, meta=meta)


@data_loader(
    "6dFGS-split", identifier=identify_6dfgs_split_fits, dtype=Spectrum1D,
    extensions=["fit", "fits"], priority=10,
)
def sixdfgs_split_fits_loader(file_obj, **kwargs):
    """
    Load the split variant of a 6dF Galaxy Survey (6dFGS) file.

    6dFGS used the Six-degree Field instrument on the UK Schmidt Telescope
    (UKST) at the Siding Spring Observatory (SSO) near Coonabarabran,
    Australia. Further details can be found at http://www.6dfgs.net/, or
    https://docs.datacentral.org.au/6dfgs/. Catalogues and spectra were
    produced, with the spectra being provided as both fits tables and as fits
    images. This loads the split variants of the spectra, which have been
    produced by either the GAMA team or Data Central.

    Parameters
    ----------
    file_obj: str, file-like or HDUList
         FITS file name, object (provided from name by Astropy I/O Registry),
         or HDUList (as resulting from astropy.io.fits.open()).

    Returns
    -------
    data: Spectrum1D
        The 6dF spectrum that is represented by the data in this file.
    """

    with read_fileobj_or_hdulist(file_obj, **kwargs) as hdulist:
        spec = _load_single_6dfgs_hdu(hdulist[0])

    return spec


@data_loader(
    "6dFGS-combined", identifier=identify_6dfgs_combined_fits,
    dtype=SpectrumList, extensions=["fit", "fits"], priority=10,
)
def sixdfgs_combined_fits_loader(file_obj, **kwargs):
    """
    Load the combined variant of a 6dF Galaxy Survey (6dFGS) file.

    6dFGS used the Six-degree Field instrument on the UK Schmidt Telescope
    (UKST) at the Siding Spring Observatory (SSO) near Coonabarabran,
    Australia. Further details can be found at http://www.6dfgs.net/, or
    https://docs.datacentral.org.au/6dfgs/. Catalogues and spectra were
    produced, with the spectra being provided as both fits tables and as fits
    images. This loads the combined variant of the spectra.

    Parameters
    ----------
    file_obj: str, file-like or HDUList
         FITS file name, object (provided from name by Astropy I/O Registry),
         or HDUList (as resulting from astropy.io.fits.open()).

    Returns
    -------
    data: SpectrumList
        The 6dF spectra that are represented by the data in this file.
    """

    with read_fileobj_or_hdulist(file_obj, **kwargs) as hdulist:
        specs = SpectrumList([_load_single_6dfgs_hdu(hdu) for hdu in hdulist[5:]])

    return specs


def _load_single_6dfgs_hdu(hdu):
    """
    Helper function to handle loading spectra from a single 6dFGS HDU
    """

    header = hdu.header
    w = WCS(naxis=1)
    w.wcs.crpix[0] = header["CRPIX1"]
    w.wcs.crval[0] = header["CRVAL1"]
    w.wcs.cdelt[0] = header["CDELT1"]
    w.wcs.cunit[0] = header["CUNIT1"]
    if w.wcs.cunit[0] == ANGSTROMS:
        w.wcs.cunit[0] = Unit("Angstrom")
    meta = {"header": header}
    flux = hdu.data[0] * Unit("count") / w.wcs.cunit[0]
    uncertainty = VarianceUncertainty(hdu.data[1])

    sky_flux = hdu.data[2] * Unit("count") / w.wcs.cunit[0]
    sky_meta = {"header": header}
    sky_spec = Spectrum1D(flux=sky_flux, wcs=w, meta=sky_meta)
    meta["sky"] = sky_spec

    return Spectrum1D(flux=flux, wcs=w, meta=meta, uncertainty=uncertainty)
