

from astropy.io import fits
import astropy.units as u
from ...spectra import Spectrum1D
from ..registers import data_loader
from astropy.nddata import InverseVariance


__all__ = ["identify_manga_cube", "identify_manga_rss", "manga_cube_loader", "manga_rss_loader"]


def _identify_sdss_fits(filename):
    """
    Check whether the given file is a SDSS data product.
    """
    try:
        with fits.open(filename, memmap=True) as hdulist:
            return hdulist[0].header["TELESCOP"] == "SDSS 2.5-M"
    except Exception:
        return False


def identify_manga_cube(origin, *args, **kwargs):
    """
    Check whether the given file is a MaNGA CUBE.
    """
    is_sdss = _identify_sdss_fits(args[0])
    with fits.open(args[0], memmap=True) as hdulist:
        return (is_sdss and "FLUX" in hdulist and hdulist[1].header['INSTRUME'] == 'MaNGA'
                and hdulist[1].header["NAXIS"] == 3)


def identify_manga_rss(origin, *args, **kwargs):
    """
    Check whether the given file is a MaNGA RSS.
    """
    is_sdss = _identify_sdss_fits(args[0])
    with fits.open(args[0], memmap=True) as hdulist:
        return (is_sdss and "FLUX" in hdulist and hdulist[1].header['INSTRUME'] == 'MaNGA'
                and hdulist[1].header["NAXIS"] == 2)


@data_loader("MaNGA cube", identifier=identify_manga_cube, dtype=Spectrum1D,
             extensions=['fits'])
def manga_cube_loader(file_obj, **kwargs):
    """
    Loader for MaNGA 3D rectified spectral data in FITS format.

    Parameters
    ----------
    filename : str
        The path to the FITS file

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
    filename : str
        The path to the FITS file

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
    unit = u.Unit('1e-17 erg / (Angstrom cm2 s)')
    if per_unit:
        unit = unit / per_unit

    hdr = hdulist['PRIMARY'].header
    wave = hdulist['WAVE'].data * u.angstrom

    if transpose:
        flux = hdulist['FLUX'].data.T * unit
        ivar = InverseVariance(hdulist["IVAR"].data.T)
    else:
        flux = hdulist['FLUX'].data * unit
        ivar = InverseVariance(hdulist["IVAR"].data)

    return Spectrum1D(flux=flux, meta={'header': hdr}, spectral_axis=wave,
                      uncertainty=ivar)
