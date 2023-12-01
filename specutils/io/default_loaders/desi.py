"""
Loader for DESI spectrum files: spectra_ and coadd_.

* spectra_ files contain all observations of an object in a given region.
* coadd_ files contain one, co-added spectrum of an object in a given region.
  The coaddition is performed across observations, but *not* across different
  cameras.

.. _spectra: https://desidatamodel.readthedocs.io/en/latest/DESI_SPECTRO_REDUX/SPECPROD/healpix/SURVEY/PROGRAM/PIXGROUP/PIXNUM/spectra-SURVEY-PROGRAM-PIXNUM.html
.. _coadd: https://desidatamodel.readthedocs.io/en/latest/DESI_SPECTRO_REDUX/SPECPROD/healpix/SURVEY/PROGRAM/PIXGROUP/PIXNUM/coadd-SURVEY-PROGRAM-PIXNUM.html
"""
import re
import warnings

from astropy.table import Table
from astropy.units import Unit
from astropy.nddata import InverseVariance
from astropy.utils.exceptions import AstropyUserWarning

import numpy as np

from ...spectra import Spectrum1D, SpectrumList
from ..registers import data_loader
from ..parsing_utils import read_fileobj_or_hdulist

__all__ = ['spectra_identify', 'coadd_identify'
           'spectra_loader', 'coadd_loader']

#
# These are reserved for future use.
#
_spectra_healpix_pattern = re.compile(r'spectra-(cmx|main|special|sv1|sv2|sv3)-(backup|bright|dark|other)-[0-9]+\.fits')
_coadd_healpix_pattern = re.compile(r'coadd-(cmx|main|special|sv1|sv2|sv3)-(backup|bright|dark|other)-[0-9]+\.fits')


def spectra_identify(origin, *args, **kwargs):
    """
    Check whether given input is FITS and has DESI stuff.
    This is used for Astropy I/O Registry.
    """
    # Test if fits has extension of type BinTable and check for spec-specific keys
    # with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
    #     return (hdulist[0].header.get('TELESCOP') == 'SDSS 2.5-M' and
    #             hdulist[0].header.get('FIBERID', 0) > 0 and
    #             len(hdulist) > 1 and
    #             (isinstance(hdulist[1], fits.BinTableHDU) and
    #              hdulist[1].header.get('TTYPE3').lower() == 'ivar'))
    return True


def coadd_identify(origin, *args, **kwargs):
    """
    Check whether given input is FITS and has DESI stuff.
    This is used for Astropy I/O Registry.
    """
    # Test if fits has extension of type BinTable and check for spec-specific keys
    # with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
    #     return (hdulist[0].header.get('TELESCOP') == 'SDSS 2.5-M' and
    #             hdulist[0].header.get('FIBERID', 0) > 0 and
    #             len(hdulist) > 1 and
    #             (isinstance(hdulist[1], fits.BinTableHDU) and
    #              hdulist[1].header.get('TTYPE3').lower() == 'ivar'))
    return True


@data_loader(
    label="DESI spectra", identifier=spectra_identify, extensions=['fits'],
    priority=10,
)
def spectra_loader(file_obj, **kwargs):
    """Loader for DESI spectra_ files.

    .. _spectra: https://desidatamodel.readthedocs.io/en/latest/DESI_SPECTRO_REDUX/SPECPROD/healpix/SURVEY/PROGRAM/PIXGROUP/PIXNUM/spectra-SURVEY-PROGRAM-PIXNUM.html

    Any keyword arguments are passed to
    :func:`~specutils.io.parsing_utils.read_fileobj_or_hdulist`.

    Parameters
    ----------
    file_obj: str, file-like, or HDUList
          FITS file name, object (provided from name by Astropy I/O Registry),
          or HDUList (as resulting from astropy.io.fits.open()).

    Returns
    -------
    SpectrumList
        Each spectrograph arm, or 'band' is represented as a Spectrum1D
        object in the SpectrumList.
    """
    return _read_desi(file_obj, **kwargs)


@data_loader(
    label="DESI coadd", identifier=coadd_identify, extensions=['fits'],
    priority=10,
)
def coadd_loader(file_obj, **kwargs):
    """Loader for DESI coadd_ files.

    .. _coadd: https://desidatamodel.readthedocs.io/en/latest/DESI_SPECTRO_REDUX/SPECPROD/healpix/SURVEY/PROGRAM/PIXGROUP/PIXNUM/coadd-SURVEY-PROGRAM-PIXNUM.html

    Any keyword arguments are passed to
    :func:`~specutils.io.parsing_utils.read_fileobj_or_hdulist`.

    Parameters
    ----------
    file_obj: str, file-like, or HDUList
          FITS file name, object (provided from name by Astropy I/O Registry),
          or HDUList (as resulting from astropy.io.fits.open()).

    Returns
    -------
    SpectrumList
        Each spectrograph arm, or 'band' is represented as a Spectrum1D
        object in the SpectrumList.
    """
    return _read_desi(file_obj, **kwargs)


def _read_desi(file_obj, **kwargs):
    """Read DESI data from a FITS file.

    This contains the common, low-level code for reading spectra and coadd files.

    Any keyword arguments are passed to
    :func:`~specutils.io.parsing_utils.read_fileobj_or_hdulist`.

    Parameters
    ----------
    file_obj: str, file-like, or HDUList
          FITS file name, object (provided from name by Astropy I/O Registry),
          or HDUList (as resulting from astropy.io.fits.open()).

    Returns
    -------
    SpectrumList
        Each spectrograph arm, or 'band' is represented as a Spectrum1D
        object in the SpectrumList.
    """
    expected_hdus = ('PRIMARY', 'FIBERMAP', 'EXP_FIBERMAP',
                     'B_WAVELENGTH', 'B_FLUX', 'B_IVAR', 'B_MASK', 'B_RESOLUTION',
                     'R_WAVELENGTH', 'R_FLUX', 'R_IVAR', 'R_MASK', 'R_RESOLUTION',
                     'Z_WAVELENGTH', 'Z_FLUX', 'Z_IVAR', 'Z_MASK', 'Z_RESOLUTION',
                     'SCORES', 'EXTRA_CATALOG')
    tables = ('FIBERMAP', 'EXP_FIBERMAP', 'SCORES', 'EXTRA_CATALOG')
    sl = SpectrumList()
    meta_zero = dict()
    band_data = dict()
    with read_fileobj_or_hdulist(file_obj, **kwargs) as hdulist:
        for k, hdu in enumerate(hdulist):
            if 'EXTNAME' in hdu.header:
                extname = hdu.header['EXTNAME']
            else:
                if k == 0:
                    extname = 'PRIMARY'
                else:
                    warnings.warn(f"HDU {k:d} has no EXTNAME and will not be read!",
                                  AstropyUserWarning)
                    extname = 'UNKNOWN'
            if extname not in expected_hdus:
                if extname != 'UNKNOWN':  # We already warned about UNKNOWN above.
                    warnings.warn(f"HDU {k:d}, EXTNAME='{extname}' is unexpected and will not be read!",
                                  AstropyUserWarning)
            elif extname == 'PRIMARY':
                meta_zero['header'] = hdu.header.copy()
            elif extname in tables:
                meta_zero[extname.lower()] = Table(hdu.data, copy=True)
            else:
                foo = extname.split('_')[0]
                band = foo[0].lower()
                datatype = foo[1]
                if band not in band_data:
                    band_data[band] = {'band': band}
                if datatype == 'WAVELENGTH':
                    band_data[band]['spectral_axis'] = hdu.data.copy() * Unit(hdu.header['BUNIT'])
                if datatype == 'FLUX':
                    meta_zero['single'] = (hdu.data.dtype == np.dtype('>f4'))
                    band_data[band]['flux'] = hdu.data.copy() * Unit(hdu.header['BUNIT'])
                if datatype == 'IVAR':
                    band_data[band]['uncertainty'] = InverseVariance(hdu.data.copy() * Unit(hdu.header['BUNIT']))
                if datatype == 'MASK':
                    band_data[band]['meta'] = {'int_mask': hdu.data.copy()}
                    band_data[band]['mask'] = band_data[band]['meta']['int_mask'] != 0
                if datatype == 'RESOLUTION':
                    pass
    meta_zero['bands'] = sorted(band_data.keys())
    for i, band in enumerate(meta_zero['bands']):
        if i == 0:
            for key, value in meta_zero.items():
                band_data[band]['meta'][key] = value
        sl.append(Spectrum1D(**(band_data[band])))
    return sl
