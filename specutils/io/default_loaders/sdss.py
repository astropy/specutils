"""
Loader for SDSS individual spectrum files: spec_ files.

.. _spec: https://data.sdss.org/datamodel/files/BOSS_SPECTRO_REDUX/RUN2D/spectra/PLATE4/spec.html
"""
import os
import re
import _io

from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.units import Unit, def_unit
from astropy.nddata import StdDevUncertainty, InverseVariance

import numpy as np

from specutils.io.registers import data_loader, custom_writer
from specutils import Spectrum1D

__all__ = ['spec_identify', 'spSpec_identify',
           'spec_loader', 'spSpec_loader']

_spSpec_pattern = re.compile(r'spSpec-\d{5}-\d{4}-\d{3}\.fit')
_spec_pattern = re.compile(r'spec-\d{4,5}-\d{5}-\d{4}\.fits')


def spec_identify(origin, *args, **kwargs):
    """
    Check whether given input is FITS and has SDSS-III/IV spec type
    BINTABLE in first extension. This is used for Astropy I/O Registry.
    """
    if isinstance(args[2], fits.hdu.hdulist.HDUList):
        hdulist = args[2]
    elif isinstance(args[2], _io.BufferedReader):
        hdulist = fits.open(args[2])
    else:
        hdulist = fits.open(args[0], **kwargs)

    # Test if fits has extension of type BinTable and check for spec-specific keys
    is_sdss = (hdulist[0].header['TELESCOP'] == 'SDSS 2.5-M' and
               hdulist[0].header.get('FIBERID', 0) > 0 and
               len(hdulist) > 1 and
               isinstance(hdulist[1], fits.BinTableHDU) and
               hdulist[1].header['TTYPE3'] == 'ivar')

    if not isinstance(args[2], (fits.hdu.hdulist.HDUList, _io.BufferedReader)):
        hdulist.close()
    return is_sdss


def spSpec_identify(origin, *args, **kwargs):
    """
    Check whether given input is FITS with SDSS-I/II spSpec tyamepe data.
    This is used for Astropy I/O Registry.
    """
    if isinstance(args[2], fits.hdu.hdulist.HDUList):
        hdulist = args[2]
    elif isinstance(args[2], _io.BufferedReader):
        hdulist = fits.open(args[2])
    else:
        hdulist = fits.open(args[0])

    # Test telescope keyword and check if primary HDU contains data
    # consistent with spSpec format
    is_sdss = (hdulist[0].header['TELESCOP'] == 'SDSS 2.5-M' and
               hdulist[0].header.get('FIBERID', 0) > 0 and
               isinstance(hdulist[0].data, np.ndarray) and
               hdulist[0].data.shape[0] == 5)

    if not isinstance(args[2], (fits.hdu.hdulist.HDUList, _io.BufferedReader)):
        hdulist.close()
    return is_sdss


def spPlate_identify(origin, *args, **kwargs):
    """
    Check whether given input is FITS with SDSS spPlate fibre spectral data.
    This is used for Astropy I/O Registry.
    """
    if isinstance(args[2], fits.hdu.hdulist.HDUList):
        hdulist = args[2]
    elif isinstance(args[2], _io.BufferedReader):
        hdulist = fits.open(args[2])
    else:
        hdulist = fits.open(args[0], **kwargs)

    # Test telescope keyword and check if primary HDU contains data
    # consistent with spSpec format
    is_sdss = (hdulist[0].header['TELESCOP'] == 'SDSS 2.5-M' and
               hdulist[0].header.get('FIBERID', 0) <= 0 and
               isinstance(hdulist[0].data, np.ndarray) and
               hdulist[0].data.shape[0] > 5)

    if not isinstance(args[2], (fits.hdu.hdulist.HDUList, _io.BufferedReader)):
        hdulist.close()
    return is_sdss


@data_loader(label="SDSS-III/IV spec", identifier=spec_identify, extensions=['fits'])
def spec_loader(file_obj, **kwargs):
    """
    Loader for SDSS-III/IV optical spectrum "spec" files.

    Parameters
    ----------
    file_obj: str, file-like, or HDUList
          FITS file name, object (provided from name by Astropy I/O Registry),
          or HDUList (as resulting from astropy.io.fits.open()).

    Returns
    -------
    data: Spectrum1D
        The spectrum that is represented by the 'loglam' (wavelength) and 'flux'
        data columns in the BINTABLE extension of the FITS `file_obj`.
    """
    if isinstance(file_obj, fits.hdu.hdulist.HDUList):
        hdulist = file_obj
    elif fits.util.fileobj_closed(file_obj):
        hdulist = fits.open(file_obj.name, **kwargs)
    else:
        hdulist = fits.open(file_obj, **kwargs)

    header = hdulist[0].header
    name = header.get('NAME')
    meta = {'header': header}

    bunit = header.get('BUNIT', '1e-17 erg / (Angstrom cm2 s)')
    if 'Ang' in bunit and 'strom' not in bunit:
        bunit = bunit.replace('Ang', 'Angstrom')
    flux_unit = Unit(bunit)

    # spectrum is in HDU 1
    flux = hdulist[1].data['flux'] * flux_unit

    uncertainty = InverseVariance(hdulist[1].data['ivar'] / flux_unit**2)

    dispersion = 10**hdulist[1].data['loglam']
    dispersion_unit = Unit('Angstrom')

    mask = hdulist[1].data['and_mask'] != 0

    if not isinstance(file_obj, fits.hdu.hdulist.HDUList):
        hdulist.close()

    return Spectrum1D(flux=flux, spectral_axis=dispersion * dispersion_unit,
                      uncertainty=uncertainty, meta=meta, mask=mask)


@data_loader(label="SDSS-I/II spSpec", identifier=spSpec_identify, extensions=['fit', 'fits'])
def spSpec_loader(file_obj, **kwargs):
    """
    Loader for SDSS-I/II spSpec files.

    Parameters
    ----------
    file_obj: str, file-like, or HDUList
           FITS file name, object (provided from name by Astropy I/O Registry),
           or HDUList (as resulting from astropy.io.fits.open()).

    Returns
    -------
    data: Spectrum1D
        The spectrum that is represented by the wavelength solution from the
        header WCS and data array of the primary HDU.
    """
    if isinstance(file_obj, fits.hdu.hdulist.HDUList):
        hdulist = file_obj
    elif fits.util.fileobj_closed(file_obj):
        hdulist = fits.open(file_obj.name, **kwargs)
    else:
        hdulist = fits.open(file_obj, **kwargs)

    header = hdulist[0].header
    # name = header.get('NAME')
    meta = {'header': header}
    wcs = WCS(header).dropaxis(1)

    bunit = header.get('BUNIT', '1e-17 erg / (Angstrom cm2 s)')
    # fix mutilated flux unit
    bunit = bunit.replace('/cm/s/Ang', '/ (Angstrom cm2 s)')
    if 'Ang' in bunit and 'strom' not in bunit:
        bunit = bunit.replace('Ang', 'Angstrom')
    flux_unit = Unit(bunit)
    flux = hdulist[0].data[0, :] * flux_unit

    uncertainty = StdDevUncertainty(hdulist[0].data[2, :] * flux_unit)

    # dispersion along NAXIS1 from the WCS
    dispersion = wcs.pixel_to_world(np.arange(flux.shape[0]))
    # convert out of logspace (default for spSpec/spPlate spectra)?
    if header.get('DC-Flag', 1) == 1:
        dispersion = 10**dispersion
    dispersion_unit = Unit('Angstrom')

    mask = hdulist[0].data[3, :] != 0

    if not isinstance(file_obj, fits.hdu.hdulist.HDUList):
        hdulist.close()

    return Spectrum1D(flux=flux, spectral_axis=dispersion * dispersion_unit,
                      uncertainty=uncertainty, meta=meta, mask=mask)


@data_loader(label="SDSS spPlate", identifier=spPlate_identify, extensions=['fits'])
def spPlate_loader(file_obj, limit=None, **kwargs):
    """
    Loader for SDSS spPlate files, reading flux spectra from all fibres into single array.

    Parameters
    ----------
    file_obj: str, file-like, or HDUList
           FITS file name, object (provided from name by Astropy I/O Registry),
           or HDUList (as resulting from astropy.io.fits.open()).

    limit : :class:`int`, optional
        If set, only return the first `limit` spectra in `flux` array.

    Returns
    -------
    Spectrum1D
        The spectra represented by the wavelength solution from the header WCS
        and the data array of the primary HDU (typically 640 along dimension 1).
    """
    if isinstance(file_obj, fits.hdu.hdulist.HDUList):
        hdulist = file_obj
    elif fits.util.fileobj_closed(file_obj):
        hdulist = fits.open(file_obj.name, **kwargs)
    else:
        hdulist = fits.open(file_obj, **kwargs)

    header = hdulist[0].header
    meta = {'header': header}
    wcs = WCS(header).dropaxis(1)
    if limit is None:
        limit = header['NAXIS2']

    bunit = header.get('BUNIT', '1e-17 erg / (Angstrom cm2 s)')
    if 'Ang' in bunit and 'strom' not in bunit:
        bunit = bunit.replace('Ang', 'Angstrom')
    flux_unit = Unit(bunit)
    flux = hdulist[0].data[0:limit, :] * flux_unit
    uncertainty = InverseVariance(hdulist[1].data[0:limit, :] / flux_unit**2)

    # dispersion along NAXIS1 from the WCS
    wcs = WCS(header).dropaxis(1)
    dispersion = wcs.pixel_to_world(np.arange(flux.shape[-1]))
    # convert out of logspace (default for spSpec/spPlate spectra)?
    if header.get('DC-Flag', 1) == 1:
        dispersion = 10**dispersion
    dispersion_unit = Unit('Angstrom')

    mask = hdulist[2].data[0:limit, :] != 0
    meta['plugmap'] = Table.read(hdulist[5])[0:limit]

    if not isinstance(file_obj, fits.hdu.hdulist.HDUList):
        hdulist.close()

    return Spectrum1D(flux=flux, spectral_axis=dispersion*dispersion_unit,
                      uncertainty=uncertainty, meta=meta, mask=mask)
