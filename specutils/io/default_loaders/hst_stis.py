from copy import deepcopy
import numpy as np
import astropy.units as u
from astropy.nddata import StdDevUncertainty

from ...spectra import Spectrum1D, SpectrumCollection
from ..registers import data_loader
from ..parsing_utils import read_fileobj_or_hdulist

__all__ = ['identify_stis_1storder_and_singleread',
           'identify_stis_multiorder_or_multiread',
           'stis_single_spectrum_loader',
           'stis_multi_spectrum_loader',]

SDQFLAGS_DEFAULT = 31743  # = 0b11110111111111
SCALAR_COLUMNS = ['SPORDER', 'NELEM', 'A2CENTER', 'EXTRSIZE', 'MAXSRCH', 'BK1SIZE',
                  'BK2SIZE', 'BK1OFFST', 'BK2OFFST', 'OFFSET',]
FLUX_UNIT = u.erg / (u.s * u.cm**2 * u.Angstrom)
DISP_UNIT = u.Angstrom


def identify_stis_1storder_and_singleread(origin, *args, **kwargs):
    """Check whether given file contains first-order HST/STIS 1D spectral data from a
    single read (ext).
    """
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        return (hdulist[0].header.get('TELESCOP', '') == 'HST') and \
               (hdulist[0].header.get('INSTRUME', '') == 'STIS') and \
               (hdulist[0].header.get('X1DCORR', '') == 'COMPLETE') and \
               (hdulist[0].header.get('NEXTEND', -1) == 1) and \
               (len(hdulist[1].data) == 1)


def identify_stis_multiorder_or_multiread(origin, *args, **kwargs):
    """Check whether given file contains HST/STIS 1D spectral data from multiple
    echelle orders and/or multiple reads (exts).
    """
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        return (hdulist[0].header.get('TELESCOP', '') == 'HST') and \
               (hdulist[0].header.get('INSTRUME', '') == 'STIS') and \
               (hdulist[0].header.get('X1DCORR', '') == 'COMPLETE') and \
               ((hdulist[0].header.get('NEXTEND', -1) > 1) or (len(hdulist[1].data) > 1))


@data_loader(
    label="HST/STIS", identifier=identify_stis_1storder_and_singleread,
    extensions=['FITS', 'FIT', 'fits', 'fit'], priority=10,
    dtype=Spectrum1D,
)
def stis_single_spectrum_loader(file_obj, **kwargs):
    """Load STIS spectral data from the MAST archive into a spectrum object.

    Parameters
    ----------
    file_obj: str, file-like, or HDUList
          FITS file name, object (provided from name by Astropy I/O Registry),
          or HDUList (as resulting from astropy.io.fits.open()).

    sdqflags: int or None
          Serious data quality flag bit-wise mask to apply.
          If None, use the default value found in the SCI ext headers.

    Returns
    -------
    data: Spectrum1D
        The spectrum that is represented by the data in this table.

    Note
    ----
    See here for the definitions of the STIS data quality (DQ) bits:
    https://hst-docs.stsci.edu/stisdhb/chapter-2-stis-data-structure/2-5-error-and-data-quality-array#id-2.5ErrorandDataQualityArray-2.5.2DataQualityFlagging
    """
    sdqflags = kwargs.pop('sdqflags', None)

    with read_fileobj_or_hdulist(file_obj, **kwargs) as hdulist:
        if len(hdulist) <= 1:
            raise ValueError('HST/STIS file is mising the SCI data extension.')
        if len(hdulist) > 2:
            raise RuntimeError('HST/STIS file is multi-extension.  '
                               'Use SpectrumCollection.read() instead.')
        if len(hdulist[1].data) > 1:
            raise RuntimeError('HST/STIS file has multiple orders (echelle data).  '
                               'Use SpectrumCollection.read() instead.')

        # Extract the single spectrum returned from the list (a Spectrum1D object):
        return _construct_Spectrum1D_list(hdulist, sdqflags=sdqflags)[0]


@data_loader(
    label="HST/STIS multi", identifier=identify_stis_multiorder_or_multiread,
    extensions=['FITS', 'FIT', 'fits', 'fit'], priority=10,
    dtype=SpectrumCollection,
)
def stis_multi_spectrum_loader(file_obj, **kwargs):
    """Load STIS spectral data from the MAST archive into a spectrum object.

    Parameters
    ----------
    file_obj: str, file-like, or HDUList
          FITS file name, object (provided from name by Astropy I/O Registry),
          or HDUList (as resulting from astropy.io.fits.open()).

    sdqflags: int or None
          Serious data quality flag bit-wise mask to apply.
          If None, use the default value found in the SCI ext headers.

    Returns
    -------
    data: SpectrumCollection
        The spectra that are represented by the data in this table.

    Note
    ----
    See here for the definitions of the STIS data quality (DQ) bits:
    https://hst-docs.stsci.edu/stisdhb/chapter-2-stis-data-structure/2-5-error-and-data-quality-array#id-2.5ErrorandDataQualityArray-2.5.2DataQualityFlagging
    """
    sdqflags = kwargs.pop('sdqflags', None)

    with read_fileobj_or_hdulist(file_obj, **kwargs) as hdulist:
        if len(hdulist) <= 1:
            raise ValueError('HST/STIS file is mising the SCI data extension.')
        if (len(hdulist) == 2) and (len(hdulist[1].data) == 1):
            raise RuntimeError('HST/STIS file is single-read and first-order.  '
                               'Use Spectrum1D.read() instead.')

        multi_spec = _construct_Spectrum1D_list(hdulist, sdqflags=sdqflags)

    return SpectrumCollection.from_spectra(multi_spec)


def _construct_Spectrum1D_list(hdulist, sdqflags=None):
    """Construct a list of Spectrum1D objects from an HST/STIS FITS HDUList.
    """
    meta = {'header_primary': hdulist[0].header}

    # Loop over all reads and orders within the dataset and accumulate Spectrum1Ds:
    multi_spec = []
    for ext in range(1, len(hdulist)):
        if sdqflags is None:
            sdqflags = hdulist[ext].header.get('SDQFLAGS', SDQFLAGS_DEFAULT)
        meta['ext'] = ext
        meta['header_sci'] = hdulist[ext].header
        for order in hdulist[ext].data:
            for scalar_column in SCALAR_COLUMNS:
                meta[scalar_column.lower()] = _nptype_to_pythontype(order[scalar_column])
            meta['sdqflags_applied'] = sdqflags
            multi_spec.append(_read_stis_order(order, meta=deepcopy(meta), sdqflags=sdqflags))

    return multi_spec


def _read_stis_order(order, meta=None, sdqflags=SDQFLAGS_DEFAULT):
    """Construct a Spectrum1D object from a single STIS order (first-order or a
    single echelle order).  Mask using SDQFLAGS.  Apply units.
    """
    return Spectrum1D(
        flux=order['FLUX'] * FLUX_UNIT,
        spectral_axis=order['WAVELENGTH'] * DISP_UNIT,
        uncertainty=StdDevUncertainty(order['ERROR'] * FLUX_UNIT),
        mask=(order['DQ'] & sdqflags) != 0,  # False where data are good
        meta=meta)


def _nptype_to_pythontype(x):
    """Truncate extra precision in decimal representation of Numpy 32-bit types when
    converting to Python-native 64-bit types.
    """
    if isinstance(x, np.float64):
        return float(x)
    elif isinstance(x, np.float32):
        return float(np.format_float_positional(x, precision=7, unique=True, trim='k'))
    elif isinstance(x, np.integer):
        return int(x)
    else:
        return float(x)
