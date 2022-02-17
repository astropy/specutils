import astropy.units as u
from astropy.nddata import InverseVariance
from astropy.wcs import WCS

from ...spectra import Spectrum1D
from ..registers import data_loader
from ..parsing_utils import read_fileobj_or_hdulist


__all__ = ["identify_manga_cube", "identify_manga_rss", "manga_cube_loader", "manga_rss_loader"]


def identify_manga_cube(origin, *args, **kwargs):
    """
    Check whether the given file is a MaNGA CUBE.
    """

    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        return (hdulist[0].header["TELESCOP"] == "SDSS 2.5-M" and "FLUX" in hdulist
                and hdulist[1].header['INSTRUME'] == 'MaNGA'
                and hdulist[1].header["NAXIS"] == 3)


def identify_manga_rss(origin, *args, **kwargs):
    """
    Check whether the given file is a MaNGA RSS.
    """

    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        return (hdulist[0].header["TELESCOP"] == "SDSS 2.5-M" and "FLUX" in hdulist
                and hdulist[1].header['INSTRUME'] == 'MaNGA'
                and hdulist[1].header["NAXIS"] == 2)


@data_loader(
    "MaNGA cube", identifier=identify_manga_cube, dtype=Spectrum1D,
    extensions=['fits'], priority=10,
)
def manga_cube_loader(file_obj, **kwargs):
    """
    Loader for MaNGA 3D rectified spectral data in FITS format.

    Parameters
    ----------
    file_obj: str, file-like, or HDUList
          FITS file name, object (provided from name by Astropy I/O Registry),
          or HDUList (as resulting from astropy.io.fits.open()).

    Returns
    -------
    Spectrum1D
        The spectrum contained in the file.
    """

    spaxel = u.Unit('spaxel', represents=u.pixel, doc='0.5" spatial pixel', parse_strict='silent')
    with read_fileobj_or_hdulist(file_obj, **kwargs) as hdulist:
        spectrum = _load_manga_spectra(hdulist, per_unit=spaxel)

    return spectrum


@data_loader(
    "MaNGA rss", identifier=identify_manga_rss, dtype=Spectrum1D,
    extensions=['fits'], priority=10,
)
def manga_rss_loader(file_obj, **kwargs):
    """
    Loader for MaNGA 2D row-stacked spectral data in FITS format.

    Parameters
    ----------
    file_obj: str, file-like, or HDUList
          FITS file name, object (provided from name by Astropy I/O Registry),
          or HDUList (as resulting from astropy.io.fits.open()).

    Returns
    -------
    Spectrum1D
        The spectrum contained in the file.
    """

    fiber = u.Unit('fiber', represents=u.pixel, doc='spectroscopic fiber', parse_strict='silent')
    with read_fileobj_or_hdulist(file_obj, **kwargs) as hdulist:
        spectrum = _load_manga_spectra(hdulist, per_unit=fiber)

    return spectrum


def _load_manga_spectra(hdulist, per_unit=None):
    """ Return a MaNGA Spectrum1D object

    Returns a Spectrum1D object for a MaNGA data files. Use the `per_unit`
    kwarg to indicate the "spaxel" or "fiber" unit for cubes and rss files,
    respectively. Note that the spectral axis will automatically be moved to
    be last during Spectrum1D initialization.

    Parameters
    ----------
    hdulist : fits.HDUList
        A MaNGA read astropy fits HDUList
    per_unit : astropy.units.Unit
        An astropy unit to divide the default flux unit by

    Returns
    -------
    Spectrum1D
        The spectrum contained in the file.
    """
    unit = u.Unit('1e-17 erg / (Angstrom cm2 s)')
    if per_unit:
        unit = unit / per_unit

    hdr = hdulist['PRIMARY'].header
    wcs = WCS(hdulist['FLUX'].header)

    flux = hdulist['FLUX'].data * unit
    ivar = InverseVariance(hdulist["IVAR"].data)
    # SDSS masks are arrays of bit values storing multiple boolean conditions.
    mask = hdulist['MASK'].data != 0

    return Spectrum1D(flux=flux, meta={'header': hdr}, wcs=wcs,
                      uncertainty=ivar, mask=mask)
