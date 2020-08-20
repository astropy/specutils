
import astropy.units as u
from astropy.io import fits
from astropy.nddata import InverseVariance

import _io
import contextlib
from ...spectra import Spectrum1D
from ..registers import data_loader

__all__ = ["identify_manga_cube", "identify_manga_rss", "manga_cube_loader", "manga_rss_loader"]


@contextlib.contextmanager
def _read_fileobj(*args, **kwargs):
    """ Context manager for reading a filename or file object

    Returns:
        an Astropy HDUList
    """
    # access the fileobj or filename arg
    # do this so identify functions are useable outside of Spectrum1d.read context
    try:
        fileobj = args[2]
    except IndexError:
        fileobj = args[0]

    if isinstance(fileobj, fits.hdu.hdulist.HDUList):
        hdulist = fileobj
    elif isinstance(fileobj, _io.BufferedReader):
        hdulist = fits.open(fileobj)
    else:
        hdulist = fits.open(fileobj, **kwargs)

    yield hdulist

    if not isinstance(fileobj, (fits.hdu.hdulist.HDUList, _io.BufferedReader)):
        hdulist.close()


def identify_manga_cube(origin, *args, **kwargs):
    """
    Check whether the given file is a MaNGA CUBE.
    """

    with _read_fileobj(*args, **kwargs) as hdulist:
        return (hdulist[0].header["TELESCOP"] == "SDSS 2.5-M" and "FLUX" in hdulist
                and hdulist[1].header['INSTRUME'] == 'MaNGA'
                and hdulist[1].header["NAXIS"] == 3)


def identify_manga_rss(origin, path, fileobj, *args, **kwargs):
    """
    Check whether the given file is a MaNGA RSS.
    """

    with _read_fileobj(*args, **kwargs) as hdulist:
        return (hdulist[0].header["TELESCOP"] == "SDSS 2.5-M" and "FLUX" in hdulist
                and hdulist[1].header['INSTRUME'] == 'MaNGA'
                and hdulist[1].header["NAXIS"] == 2)


@data_loader("MaNGA cube", identifier=identify_manga_cube, dtype=Spectrum1D,
             extensions=['fits'])
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

    if isinstance(file_obj, fits.hdu.hdulist.HDUList):
        hdulist = file_obj
    else:
        hdulist = fits.open(file_obj, **kwargs)

    spaxel = u.Unit('spaxel', represents=u.pixel, doc='0.5" spatial pixel', parse_strict='silent')
    spectrum = _load_manga_spectra(hdulist, per_unit=spaxel, transpose=True)

    if not isinstance(file_obj, fits.hdu.hdulist.HDUList):
        hdulist.close()

    return spectrum


@data_loader("MaNGA rss", identifier=identify_manga_rss, dtype=Spectrum1D,
             extensions=['fits'])
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

    if isinstance(file_obj, fits.hdu.hdulist.HDUList):
        hdulist = file_obj
    else:
        hdulist = fits.open(file_obj, **kwargs)

    fiber = u.Unit('fiber', represents=u.pixel, doc='spectroscopic fiber', parse_strict='silent')
    spectrum = _load_manga_spectra(hdulist, per_unit=fiber)

    if not isinstance(file_obj, fits.hdu.hdulist.HDUList):
        hdulist.close()

    return spectrum


def _load_manga_spectra(hdulist, per_unit=None, transpose=None):
    """ Return a MaNGA Spectrum1D object

    Returns a Spectrum1D object for a MaNGA data files.  Set
    `transpose` kwarg to True for MaNGA cubes, as they are flipped relative
    to what Spectrum1D expects.  Use the `per_unit` kwarg to indicate the
    "spaxel" or "fiber" unit for cubes and rss files, respectively.

    Parameters
    ----------
    hdulist : fits.HDUList
        A MaNGA read astropy fits HDUList
    per_unit : astropy.units.Unit
        An astropy unit to divide the default flux unit by
    transpose : bool
        If True, transpose the data arrays

    Returns
    -------
    Spectrum1D
        The spectrum contained in the file.
    """
    unit = u.Unit('1e-17 erg / (Angstrom cm2 s)')
    if per_unit:
        unit = unit / per_unit

    hdr = hdulist['PRIMARY'].header
    wave = hdulist['WAVE'].data * u.angstrom

    if transpose:
        flux = hdulist['FLUX'].data.T * unit
        ivar = InverseVariance(hdulist["IVAR"].data.T)
        # SDSS masks are arrays of bit values storing multiple boolean conditions.
        # Setting non-zero bit values to True to map to specutils standard
        mask = hdulist['MASK'].data.T != 0
    else:
        flux = hdulist['FLUX'].data * unit
        ivar = InverseVariance(hdulist["IVAR"].data)
        mask = hdulist['MASK'].data != 0

    return Spectrum1D(flux=flux, meta={'header': hdr}, spectral_axis=wave,
                      uncertainty=ivar, mask=mask)
