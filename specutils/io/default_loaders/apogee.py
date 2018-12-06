"""
Loader for APOGEE spectrum files: apVisit_, apStar_, aspcapStar_ files.

.. _apVisit: https://data.sdss.org/datamodel/files/APOGEE_REDUX/APRED_VERS/TELESCOPE/PLATE_ID/MJD5/apVisit.html
.. _apStar: https://data.sdss.org/datamodel/files/APOGEE_REDUX/APRED_VERS/APSTAR_VERS/TELESCOPE/LOCATION_ID/apStar.html
.. _aspcapStar: https://data.sdss.org/datamodel/files/APOGEE_REDUX/APRED_VERS/APSTAR_VERS/ASPCAP_VERS/RESULTS_VERS/LOCATION_ID/aspcapStar.html
"""
import os

from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.units import Unit, def_unit
from astropy.nddata import StdDevUncertainty

import numpy as np

from specutils.io.registers import data_loader, custom_writer
from specutils import Spectrum1D

__all__ = ['apVisit_identify', 'apStar_identify', 'aspcapStar_identify',
           'apVisit_loader', 'apStar_loader', 'aspcapStar_loader']


def apVisit_identify(origin, *args, **kwargs):
    """
    Check whether given filename is FITS. This is used for Astropy I/O
    Registry.
    """
    return (isinstance(args[0], str) and
            args[0].lower().split('.')[-1] == 'fits' and
            args[0].startswith('apVisit'))


def apStar_identify(origin, *args, **kwargs):
    """
    Check whether given filename is FITS. This is used for Astropy I/O
    Registry.
    """
    return (isinstance(args[0], str) and
            args[0].lower().split('.')[-1] == 'fits' and
            args[0].startswith('apStar'))


def aspcapStar_identify(origin, *args, **kwargs):
    """
    Check whether given filename is FITS. This is used for Astropy I/O
    Registry.
    """
    return (isinstance(args[0], str) and
            args[0].lower().split('.')[-1] == 'fits' and
            args[0].startswith('aspcapStar'))


@data_loader(label="APOGEE apVisit", identifier=apVisit_identify, extensions=['fits'])
def apVisit_loader(file_name, **kwargs):
    """
    Loader for APOGEE apVisit files.

    Parameters
    ----------
    file_name: str
        The path to the FITS file

    Returns
    -------
    data: Spectrum1D
        The spectrum that is represented by the data in this table.
    """
    name = os.path.basename(file_name.rstrip(os.sep)).rsplit('.', 1)[0]
    hdulist = fits.open(file_name, **kwargs)

    header = hdulist[0].header
    meta = {'header': header}

    # spectrum is stored in three rows (for three chips)
    data = np.concatenate([hdulist[1].data[0, :],
                           hdulist[1].data[1, :],
                           hdulist[1].data[2, :]])
    unit = Unit('1e-17 erg / (Angstrom cm2 s)')

    stdev = np.concatenate([hdulist[2].data[0, :],
                           hdulist[2].data[1, :],
                           hdulist[2].data[2, :]])
    uncertainty = StdDevUncertainty(stdev * unit)

    # Dispersion is not a simple function in these files.  There's a
    # look-up table instead.
    dispersion = np.concatenate([hdulist[4].data[0, :],
                                 hdulist[4].data[1, :],
                                 hdulist[4].data[2, :]])
    dispersion_unit = Unit('Angstrom')
    hdulist.close()

    return Spectrum1D(data=data * unit,
                      uncertainty=uncertainty,
                      dispersion=dispersion * dispersion_unit,
                      meta=meta)


@data_loader(label="APOGEE apStar", identifier=apStar_identify, extensions=['fits'])
def apStar_loader(file_name, **kwargs):
    """
    Loader for APOGEE apStar files.

    Parameters
    ----------
    file_name: str
        The path to the FITS file

    Returns
    -------
    data: Spectrum1D
        The spectrum that is represented by the data in this table.
    """
    name = os.path.basename(file_name.rstrip(os.sep)).rsplit('.', 1)[0]
    hdulist = fits.open(file_name, **kwargs)

    header = hdulist[0].header
    meta = {'header': header}
    wcs = WCS(hdulist[1].header)

    data = hdulist[1].data[0, :]  # spectrum in the first row of the first extension
    unit = Unit('1e-17 erg / (Angstrom cm2 s)')

    uncertainty = StdDevUncertainty(hdulist[2].data[0, :])

    # dispersion from the WCS but convert out of logspace
    # dispersion = 10**wcs.all_pix2world(np.arange(data.shape[0]), 0)[0]
    dispersion = 10**wcs.all_pix2world(np.vstack((np.arange(data.shape[0]),
                                                  np.zeros((data.shape[0],)))).T,
                                       0)[:, 0]
    dispersion_unit = Unit('Angstrom')
    hdulist.close()

    return Spectrum1D(data=data * unit,
                      uncertainty=uncertainty,
                      dispersion=dispersion * dispersion_unit,
                      meta=meta,
                      wcs=wcs)


@data_loader(label="APOGEE aspcapStar", identifier=aspcapStar_identify, extensions=['fits'])
def aspcapStar_loader(file_name, **kwargs):
    """
    Loader for APOGEE aspcapStar files.

    Parameters
    ----------
    file_name: str
        The path to the FITS file

    Returns
    -------
    data: Spectrum1D
        The spectrum that is represented by the data in this table.
    """
    name = os.path.basename(file_name.rstrip(os.sep)).rsplit('.', 1)[0]
    hdulist = fits.open(file_name, **kwargs)

    header = hdulist[0].header
    meta = {'header': header}
    wcs = WCS(hdulist[1].header)

    data = hdulist[1].data # spectrum in the first extension
    unit = def_unit('arbitrary units')

    uncertainty = StdDevUncertainty(hdulist[2].data)

    # dispersion from the WCS but convert out of logspace
    dispersion = 10**wcs.all_pix2world(np.arange(data.shape[0]), 0)[0]
    dispersion_unit = Unit('Angstrom')
    hdulist.close()

    return Spectrum1D(data=data * unit,
                      uncertainty=uncertainty,
                      dispersion=dispersion * dispersion_unit,
                      meta=meta,
                      wcs=wcs)
