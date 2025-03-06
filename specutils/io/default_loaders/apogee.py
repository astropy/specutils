"""
Loader for APOGEE spectrum files: apVisit_, apStar_, aspcapStar_ files.

.. _apVisit: https://data.sdss.org/datamodel/files/APOGEE_REDUX/APRED_VERS/TELESCOPE/PLATE_ID/MJD5/apVisit.html
.. _apStar: https://data.sdss.org/datamodel/files/APOGEE_REDUX/APRED_VERS/APSTAR_VERS/TELESCOPE/LOCATION_ID/apStar.html
.. _aspcapStar: https://data.sdss.org/datamodel/files/APOGEE_REDUX/APRED_VERS/APSTAR_VERS/ASPCAP_VERS/RESULTS_VERS/LOCATION_ID/aspcapStar.html
"""
import numpy as np
from astropy.nddata import StdDevUncertainty
from astropy.units import Unit, def_unit
from astropy.wcs import WCS

from ...spectra import Spectrum1D
from ..parsing_utils import read_fileobj_or_hdulist
from ..registers import data_loader

__all__ = ['apVisit_identify', 'apStar_identify', 'aspcapStar_identify',
           'apVisit_loader', 'apStar_loader', 'aspcapStar_loader']


def apVisit_identify(origin, *args, **kwargs):
    """
    Check whether given file is FITS. This is used for Astropy I/O
    Registry.
    """
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        # Test if fits has extension of type BinTable and check for
        # apVisit-specific keys
        return (
            "apogee" in hdulist[0].header.get("SURVEY", "none")
            and "APOGEE" in "".join(hdulist[0].header["HISTORY"])
            and len(hdulist) > 4
            and hdulist[1].header.get("BUNIT", "none").startswith("Flux")
            and hdulist[2].header.get("BUNIT", "none").startswith("Flux")
            and hdulist[4].header.get("BUNIT", "none").startswith("Wavelen")
        )


def apStar_identify(origin, *args, **kwargs):
    """
    Check whether given file is FITS. This is used for Astropy I/O
    Registry.
    """
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        # Test if fits has extension of type BinTable and check for
        # apogee-specific keys + keys for individual apVisits
        return "APOGEE" in "".join(hdulist[0].header["HISTORY"]) and hdulist[0].header.get("SFILE1", "none").startswith(
            ("apVisit", "asVisit")
        )


def aspcapStar_identify(origin, *args, **kwargs):
    """
    Check whether given file is FITS. This is used for Astropy I/O
    Registry.
    """
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        # Test if fits has extension of type BinTable and check for
        # aspcapStar-specific keys
        return (
            len(hdulist) > 4
            and hdulist[1].header.get("NAXIS1", 0) > 8000
            and hdulist[2].header.get("NAXIS1", 0) > 8000
            and "ASPCAPFLAG" in hdulist[-1].header.values()
        )


@data_loader(
    label="APOGEE apVisit", identifier=apVisit_identify, extensions=['fits'],
    priority=10,
)
def apVisit_loader(file_obj, **kwargs):
    """
    Loader for APOGEE apVisit files.

    Parameters
    ----------
    file_obj: str, file-like or HDUList
        FITS file name, object (provided from name by Astropy I/O Registry),
        or HDUList (as resulting from astropy.io.fits.open()).

    Returns
    -------
    data: Spectrum1D
        The spectrum that is represented by the data in this table.
    """

    with read_fileobj_or_hdulist(file_obj, **kwargs) as hdulist:
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

    return Spectrum1D(data=data * unit,
                      uncertainty=uncertainty,
                      spectral_axis=dispersion * dispersion_unit,
                      meta=meta)


@data_loader(
    label="APOGEE apStar", identifier=apStar_identify, extensions=['fits'],
    priority=10,
)
def apStar_loader(file_obj, **kwargs):
    """
    Loader for APOGEE apStar files.

    Parameters
    ----------
    file_obj: str, file-like, or HDUList
        FITS file name or object (provided from name by Astropy I/O Registry),
        or HDUList (as resulting from astropy.io.fits.open()).

    Returns
    -------
    data: Spectrum1D
        The spectrum that is represented by the data in this table.
    """

    with read_fileobj_or_hdulist(file_obj, **kwargs) as hdulist:
        header = hdulist[0].header
        meta = {'header': header}
        wcs = WCS(hdulist[1].header)

        data = hdulist[1].data[0, :]  # spectrum in the first row of the first extension
        unit = Unit('1e-17 erg / (Angstrom cm2 s)')

        uncertainty = StdDevUncertainty(hdulist[2].data[0, :])

    # dispersion from the WCS but convert out of logspace
    # dispersion = 10**wcs.all_pix2world(np.arange(data.shape[0]), 0)[0]
    dispersion = 10**wcs.all_pix2world(np.vstack((np.arange(data.shape[0]),
                                                  np.zeros((data.shape[0],)))).T, 0)[:, 0]
    dispersion_unit = Unit('Angstrom')

    return Spectrum1D(data=data * unit,
                      uncertainty=uncertainty,
                      spectral_axis=dispersion * dispersion_unit,
                      meta=meta)


@data_loader(
    label="APOGEE aspcapStar", identifier=aspcapStar_identify,
    extensions=['fits'], priority=10,
)
def aspcapStar_loader(file_obj, **kwargs):
    """
    Loader for APOGEE aspcapStar files.

    Parameters
    ----------
    file_obj: str or file-like
        FITS file name or object (provided from name by Astropy I/O Registry).

    Returns
    -------
    data: Spectrum1D
        The spectrum that is represented by the data in this table.
    """

    with read_fileobj_or_hdulist(file_obj, **kwargs) as hdulist:
        header = hdulist[0].header
        meta = {'header': header}
        wcs = WCS(hdulist[1].header)

        data = hdulist[1].data  # spectrum in the first extension
        unit = def_unit('arbitrary units')

        uncertainty = StdDevUncertainty(hdulist[2].data)

    # dispersion from the WCS but convert out of logspace
    dispersion = 10**wcs.all_pix2world(np.arange(data.shape[0]), 0)[0]
    dispersion_unit = Unit('Angstrom')

    return Spectrum1D(data=data * unit,
                      uncertainty=uncertainty,
                      spectral_axis=dispersion * dispersion_unit,
                      meta=meta,
                      wcs=wcs)
