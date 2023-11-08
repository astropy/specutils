"""
Loader for SDSS spectrum files: spec_, spSpec_, and spPlate_.

.. _spec: https://data.sdss.org/datamodel/files/BOSS_SPECTRO_REDUX/RUN2D/spectra/PLATE4/spec.html
.. _spSpec: http://classic.sdss.org/dr7/dm/flatFiles/spSpec.php
.. _spPlate: https://data.sdss.org/datamodel/files/BOSS_SPECTRO_REDUX/RUN2D/PLATE4/spPlate.html
"""
import re

from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.units import Unit
from astropy.nddata import StdDevUncertainty, InverseVariance

import numpy as np

from ...spectra import Spectrum1D
from ..registers import data_loader
from ..parsing_utils import read_fileobj_or_hdulist

__all__ = ['spec_identify', 'spSpec_identify', 'spPlate_identify',
           'spec_loader', 'spSpec_loader', 'spPlate_loader']

#
# These are reserved for future use.
#
_spSpec_pattern = re.compile(r'spSpec-\d{5}-\d{4}-\d{3}\.fit')
_spec_pattern = re.compile(r'spec-\d{4,5}-\d{5}-\d{4}\.fits')
_spPlate_pattern = re.compile(r'spPlate-\d{4,5}-\d{5}\.fits')


def _sdss_wcs_to_log_wcs(old_wcs):
    """
    The WCS in the SDSS files does not appear to follow the WCS standard - it
    claims to be linear, but is logarithmic in base-10.
    The wavelength is given by:
    λ = 10^(w0 + w1 * i)
    with i being the pixel index starting from 0. This formula is documented at
    https://classic.sdss.org/dr7/products/spectra/read_spSpec.php and appears to
    be the same across all SDSS-I and SDSS-II data releases (replace dr7 with
    dr<x> in the above URL).

    The FITS standard uses a natural log with a sightly different formulation,
    see WCS Paper 3 (which discusses spectral WCS).

    This function does the conversion from the SDSS WCS to FITS WCS.
    """
    w0 = old_wcs.wcs.crval[0]
    w1 = old_wcs.wcs.cd[0, 0]
    crval = 10**w0
    cdelt = crval * w1 * np.log(10)
    cunit = Unit("Angstrom")
    ctype = "WAVE-LOG"

    w = WCS(naxis=1)
    w.wcs.crval[0] = crval
    w.wcs.cdelt[0] = cdelt
    w.wcs.ctype[0] = ctype
    w.wcs.cunit[0] = cunit
    w.wcs.set()

    return w


def spec_identify(origin, *args, **kwargs):
    """
    Check whether given input is FITS and has SDSS-III/IV spec type
    BINTABLE in first extension. This is used for Astropy I/O Registry.
    """
    # Test if fits has extension of type BinTable and check for spec-specific keys
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        return (hdulist[0].header.get("TELESCOP") == "SDSS 2.5-M"
                and hdulist[0].header.get("FIBERID", 0) > 0
                and len(hdulist) > 1
                and (isinstance(hdulist[1], fits.BinTableHDU)
                     and hdulist[1].header.get("TTYPE3").lower() == "ivar"))


def spSpec_identify(origin, *args, **kwargs):
    """
    Check whether given input is FITS with SDSS-I/II spSpec tyamepe data.
    This is used for Astropy I/O Registry.
    """
    # Test telescope keyword and check if primary HDU contains data
    # consistent with spSpec format
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        return (hdulist[0].header.get("TELESCOP") == "SDSS 2.5-M"
                and hdulist[0].header.get("FIBERID", 0) > 0
                and isinstance(hdulist[0].data, np.ndarray)
                and hdulist[0].data.shape[0] == 5)


def spPlate_identify(origin, *args, **kwargs):
    """
    Check whether given input is FITS with SDSS spPlate fibre spectral data.
    This is used for Astropy I/O Registry.
    """
    # Test telescope keyword and check if primary HDU contains data
    # consistent with spSpec format
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        return (hdulist[0].header.get("TELESCOP") == "SDSS 2.5-M"
                and hdulist[0].header.get("FIBERID", 0) <= 0
                and isinstance(hdulist[0].data, np.ndarray)
                and hdulist[0].data.shape[0] > 5)


@data_loader(
    label="SDSS-III/IV spec",
    identifier=spec_identify,
    extensions=["fits"],
    priority=10,
)
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
    with read_fileobj_or_hdulist(file_obj, **kwargs) as hdulist:
        header = hdulist[0].header
        meta = {"header": header}

        bunit = header.get("BUNIT", "1e-17 erg / (Angstrom cm2 s)")
        if "Ang" in bunit and "strom" not in bunit:
            bunit = bunit.replace("Ang", "Angstrom")
        flux_unit = Unit(bunit)

        # spectrum is in HDU 1
        flux = hdulist[1].data["flux"] * flux_unit

        uncertainty = InverseVariance(hdulist[1].data["ivar"] / flux_unit**2)

        dispersion = 10**hdulist[1].data["loglam"]
        dispersion_unit = Unit("Angstrom")

        mask = hdulist[1].data["and_mask"] != 0

    return Spectrum1D(
        flux=flux,
        spectral_axis=dispersion * dispersion_unit,
        uncertainty=uncertainty,
        meta=meta,
        mask=mask,
    )


@data_loader(
    label="SDSS-I/II spSpec",
    identifier=spSpec_identify,
    extensions=["fit", "fits"],
    priority=10,
)
def spSpec_loader(file_obj, **kwargs):
    """
    Loader for SDSS-I/II spSpec files.

    The content of these files is documented at
    https://classic.sdss.org/dr7/dm/flatFiles/spSpec.php, with instructions for
    reading them at
    https://classic.sdss.org/dr7/products/spectra/read_spSpec.php.

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
    with read_fileobj_or_hdulist(file_obj, **kwargs) as hdulist:
        header = hdulist[0].header
        # name = header.get('NAME')
        meta = {"header": header}
        wcs = WCS(header).dropaxis(1)

        bunit = header.get("BUNIT", "1e-17 erg / (Angstrom cm2 s)")
        # fix mutilated flux unit
        bunit = bunit.replace("/cm/s/Ang", "/ (Angstrom cm2 s)")
        if "Ang" in bunit and "strom" not in bunit:
            bunit = bunit.replace("Ang", "Angstrom")
        flux_unit = Unit(bunit)
        flux = hdulist[0].data[0, :] * flux_unit

        uncertainty = StdDevUncertainty(hdulist[0].data[2, :] * flux_unit)

        # Fix the WCS if it is claimed to be linear
        if header.get("DC-Flag", 1) == 1:
            fixed_wcs = _sdss_wcs_to_log_wcs(wcs)
        else:
            fixed_wcs = wcs

        mask = hdulist[0].data[3, :] != 0

    return Spectrum1D(flux=flux,
                      wcs=fixed_wcs,
                      uncertainty=uncertainty,
                      meta=meta,
                      mask=mask)


@data_loader(label="SDSS spPlate",
             identifier=spPlate_identify,
             extensions=["fits"])
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
    with read_fileobj_or_hdulist(file_obj, **kwargs) as hdulist:
        header = hdulist[0].header
        meta = {"header": header}
        wcs = WCS(header).dropaxis(1)
        if limit is None:
            limit = header["NAXIS2"]

        bunit = header.get("BUNIT", "1e-17 erg / (Angstrom cm2 s)")
        if "Ang" in bunit and "strom" not in bunit:
            bunit = bunit.replace("Ang", "Angstrom")
        flux_unit = Unit(bunit)
        flux = hdulist[0].data[0:limit, :] * flux_unit
        uncertainty = InverseVariance(hdulist[1].data[0:limit, :] /
                                      flux_unit**2)

        # Fix the WCS if it is claimed to be linear
        if header.get("DC-Flag", 1) == 1:
            fixed_wcs = _sdss_wcs_to_log_wcs(wcs)
        else:
            fixed_wcs = wcs

        mask = hdulist[2].data[0:limit, :] != 0
        meta["plugmap"] = Table.read(hdulist[5])[0:limit]

    return Spectrum1D(flux=flux,
                      wcs=fixed_wcs,
                      uncertainty=uncertainty,
                      meta=meta,
                      mask=mask)
