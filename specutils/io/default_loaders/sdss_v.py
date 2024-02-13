"""Register reader functions for various spectral formats."""
from typing import Optional

import numpy as np
from astropy.units import Unit, Quantity, Angstrom
from astropy.nddata import StdDevUncertainty, InverseVariance
from astropy.io.fits import HDUList, BinTableHDU, ImageHDU

from ...spectra import Spectrum1D, SpectrumList
from ..registers import data_loader
from ..parsing_utils import read_fileobj_or_hdulist

__all__ = [
    "load_sdss_apVisit_1D",
    "load_sdss_apVisit_list",
    "load_sdss_apStar_1D",
    "load_sdss_apStar_list",
    "load_sdss_spec_1D",
    "load_sdss_spec_list",
    "load_sdss_mwm_1d",
    "load_sdss_mwm_list",
]


def apVisit_identify(origin, *args, **kwargs):
    """
    Check whether given file is FITS. This is used for Astropy I/O
    Registry.

    Updated for SDSS-V outputs.
    """
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        # Test if fits has extension of type BinTable and check for
        # apVisit-specific keys
        return (hdulist[0].header.get("SURVEY").lower() == "sdss-v"
                and len(hdulist) > 4
                and hdulist[1].header.get("BUNIT", "none").startswith("Flux")
                and hdulist[2].header.get("BUNIT", "none").startswith("Flux")
                and hdulist[4].header.get("BUNIT", "none").startswith("Wave"))


def apStar_identify(origin, *args, **kwargs):
    """
    Check whether given file is FITS. This is used for Astropy I/O
    Registry.

    Updated for SDSS-V outputs.
    """
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        # Test if fits has extension of type BinTable and check for
        # apogee-specific keys + keys for individual apVisits
        return ("APRED" in hdulist[0].header.keys() and len(hdulist) > 9
                and "NVISITS" in hdulist[0].header.keys()
                and isinstance(hdulist[1], ImageHDU)
                and hdulist[1].header.get("BUNIT", "none").startswith("Flux")
                and hdulist[2].header.get("BUNIT", "none").startswith("Err"))


def spec_sdss5_identify(origin, *args, **kwargs):
    """
    Check whether given input is FITS and has SDSS-V spec type
    BINTABLE in first extension. This is used for Astropy I/O Registry.

    Derived from SDSS-III/IV method in sdss.py. Single change (FIBERID key).
    """
    # Test if fits has extension of type BinTable and check for spec-specific keys
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        return (
            (hdulist[0].header.get("TELESCOP").lower() == "sdss 2.5-m")  # and
            # hdulist[0].header.get("OBSERVAT").lower() in ["apo", "lco"]
            and (hdulist[1].header.get("TTYPE1").lower() == "flux") and
            (hdulist[1].header.get("TTYPE2").lower() == "loglam") and
            (len(hdulist) > 1) and (isinstance(hdulist[1], BinTableHDU)) and
            (hdulist[1].header.get("TTYPE3").lower() == "ivar"))


def mwm_identify(origin, *args, **kwargs):
    """
    Check whether given input is FITS and has SDSS-V mwm type
    BINTABLE in all 4 extensions. This is used for Astropy I/O Registry.
    """
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        return (("V_ASTRA" in hdulist[0].header.keys()) and len(hdulist) > 0
                and ("SDSS_ID" in hdulist[0].header.keys())
                and (isinstance(hdulist[i], BinTableHDU) for i in range(1, 5)))


def _wcs_log_linear(naxis, cdelt, crval):
    """
    Convert WCS from log to linear.
    """
    return 10**(np.arange(naxis) * cdelt + crval)


def _fetch_flux_unit(hdu):
    """
    Fetch the flux unit by accessing the BUNIT key of the given HDU.

    Parameters
    ----------
    hdu : HDUImage
        Flux or e_flux HDUImage object.

    Returns
    -------
    flux_unit : Unit
        Flux unit as Astropy Unit object.

    """
    flux_unit = (hdu.header.get("BUNIT").replace("(", "").replace(
        ")", "").split(" ")[-2:])
    flux_unit = "".join(flux_unit)
    if "Ang" in flux_unit and "strom" not in flux_unit:
        flux_unit = flux_unit.replace("Ang", "Angstrom")

    flux_unit = flux_unit.replace("s/cm^2/Angstrom", "(s cm2 Angstrom)")

    return Unit(flux_unit)


# APOGEE files
@data_loader(
    "SDSS-V apStar",
    identifier=apStar_identify,
    dtype=Spectrum1D,
    priority=10,
    extensions=["fits"],
)
def load_sdss_apStar_1D(file_obj, idx: int = 0, **kwargs):
    """
    Load an apStar file as a Spectrum1D.

    Parameters
    ----------
    file_obj : str, file-like, or HDUList
        FITS file name, file object, or HDUList.
    idx : int
        The specified visit to load. Defaults to coadd (index 0).

    Returns
    -------
    Spectrum1D
        The spectrum contained in the file

    """
    # intialize slicer

    with read_fileobj_or_hdulist(file_obj, memmap=False, **kwargs) as hdulist:
        # Obtain spectral axis from header (HDU[1])
        wavelength = _wcs_log_linear(
            hdulist[1].header.get("NAXIS1"),
            hdulist[1].header.get("CDELT1"),
            hdulist[1].header.get("CRVAL1"),
        )
        spectral_axis = Quantity(wavelength, unit=Angstrom)

        # Obtain flux and e_flux from data (HDU[1] & HDU[2])
        flux_unit = _fetch_flux_unit(hdulist[1])
        flux = Quantity(hdulist[1].data[idx], unit=flux_unit)
        e_flux = StdDevUncertainty(hdulist[2].data[idx])

        # reduce flux array if 1D in 2D np array
        # NOTE: bypasses jdaviz specviz NotImplementedError, but could be the expected output for specutils
        if flux.shape[0] == 1:
            flux = flux[0]
            e_flux = e_flux[0]

        # Add header + bitmask
        meta = dict()
        meta["header"] = hdulist[0].header
        # SDSS masks are arrays of bit values storing multiple bool conds
        mask = hdulist[3].data[idx]
        # NOTE: specutils considers 0/False as valid values, simlar to numpy convention
        mask = mask != 0

    return Spectrum1D(
        spectral_axis=spectral_axis,
        flux=flux,
        uncertainty=e_flux,
        mask=mask,
        meta=meta,
    )


@data_loader(
    "SDSS-V apStar multi",
    identifier=apStar_identify,
    dtype=SpectrumList,
    priority=10,
    extensions=["fits"],
)
def load_sdss_apStar_list(file_obj, **kwargs):
    """
    Load an apStar file as a SpectrumList.

    Parameters
    ----------
    file_obj : str, file-like, or HDUList
        FITS file name, file object, or HDUList.

    Returns
    -------
    SpectrumList
        The spectra contained in the file
    """
    with read_fileobj_or_hdulist(file_obj, **kwargs) as hdulist:
        nvisits = hdulist[0].header.get("NVISITS")
        if nvisits <= 1:
            raise ValueError(
                "Only 1 visit in this file. Use Spectrum1D.read() instead.")
        return SpectrumList([
            load_sdss_apStar_1D(file_obj, idx=i, **kwargs)
            for i in range(nvisits)
        ])


@data_loader(
    "SDSS-V apVisit",
    identifier=apVisit_identify,
    dtype=Spectrum1D,
    priority=10,
    extensions=["fits"],
)
def load_sdss_apVisit_1D(file_obj, **kwargs):
    """
    Load an apVisit file as a Spectrum1D.

    Parameters
    ----------
    file_obj : str, file-like, or HDUList
        FITS file name, file object, or HDUList.

    Returns
    -------
    Spectrum1D
        The spectrum contained in the file
    """
    flux_unit = Unit("1e-17 erg / (Angstrom cm2 s)")

    with read_fileobj_or_hdulist(file_obj, memmap=False, **kwargs) as hdulist:
        # Fetch and flatten spectrum parameters
        spectral_axis = Quantity(hdulist[4].data.flatten(), unit=Angstrom)
        flux_unit = _fetch_flux_unit(hdulist[1])
        flux = Quantity(hdulist[1].data.flatten(), unit=flux_unit)
        e_flux = StdDevUncertainty(hdulist[2].data.flatten())

        # Get metadata and attach bitmask and MJD
        meta = dict()
        meta["header"] = hdulist[0].header
        mask = hdulist[3].data.flatten()
        # NOTE: specutils considers 0/False as valid values, simlar to numpy convention
        mask = mask != 0

    return Spectrum1D(spectral_axis=spectral_axis,
                      flux=flux,
                      mask=mask,
                      uncertainty=e_flux,
                      meta=meta)


@data_loader(
    "SDSS-V apVisit multi",
    identifier=apVisit_identify,
    dtype=SpectrumList,
    priority=10,
    extensions=["fits"],
)
def load_sdss_apVisit_list(file_obj, **kwargs):
    """
    Load an apVisit file as a SpectrumList.

    Parameters
    ----------
    file_obj : str, file-like, or HDUList
        FITS file name, file object, or HDUList.

    Returns
    -------
    SpectrumList
        The spectra from each chip contained in the file
    """
    spectra = SpectrumList()
    with read_fileobj_or_hdulist(file_obj, memmap=False, **kwargs) as hdulist:
        # Get metadata
        flux_unit = _fetch_flux_unit(hdulist[1])
        common_meta = dict()
        common_meta["header"] = hdulist[0].header

        for chip in range(hdulist[1].data.shape[0]):
            # Fetch spectral axis and flux, and E_flux
            spectral_axis = Quantity(hdulist[4].data[chip], unit=Angstrom)
            flux = Quantity(hdulist[1].data[chip], unit=flux_unit)
            e_flux = StdDevUncertainty(hdulist[2].data[chip])

            # Copy metadata for each, adding chip bitmask and MJD to each
            meta = common_meta.copy()
            mask = hdulist[3].data[chip]
            # NOTE: specutils considers 0/False as valid values, simlar to numpy convention
            mask = mask != 0

            spectra.append(
                Spectrum1D(
                    spectral_axis=spectral_axis,
                    flux=flux,
                    mask=mask,
                    uncertainty=e_flux,
                    meta=meta,
                ))

    return spectra


# BOSS REDUX products (specLite, specFull, custom coadd files, etc)

@data_loader(
    "SDSS-V spec",
    identifier=spec_sdss5_identify,
    dtype=Spectrum1D,
    priority=5,
    extensions=["fits"],
)
def load_sdss_spec_1D(file_obj, *args, hdu: Optional[int] = None, **kwargs):
    """
    Load a given BOSS spec file as a Spectrum1D object.

    Parameters
    ----------
    file_obj : str, file-like, or HDUList
        FITS file name, file object, or HDUList..
    hdu : int
        The specified HDU to load a given spectra from.

    Returns
    -------
    Spectrum1D
        The spectrum contained in the file at the specified HDU.
    """
    if hdu is None:
        # TODO: how should we handle this -- multiple things in file, but the user cannot choose.
        hdu = 1  # defaulting to coadd
        # raise ValueError("HDU not specified! Please specify a HDU to load.")
    elif hdu in [2, 3, 4]:
        raise ValueError("Invalid HDU! HDU{} is not spectra.".format(hdu))
    with read_fileobj_or_hdulist(file_obj, memmap=False, **kwargs) as hdulist:
        # directly load the coadd at HDU1
        return _load_BOSS_HDU(hdulist, hdu, **kwargs)


@data_loader(
    "SDSS-V spec multi",
    identifier=spec_sdss5_identify,
    dtype=SpectrumList,
    priority=5,
    extensions=["fits"],
)
def load_sdss_spec_list(file_obj, **kwargs):
    """
    Load a given BOSS spec file as a SpectrumList object.

    Parameters
    ----------
    file_obj : str, file-like, or HDUList
        FITS file name, file object, or HDUList..

    Returns
    -------
    SpectrumList
        The spectra contained in the file.
    """
    with read_fileobj_or_hdulist(file_obj, memmap=False, **kwargs) as hdulist:
        spectra = list()
        for hdu in range(1, len(hdulist)):
            if hdulist[hdu].name in ["SPALL", "ZALL", "ZLINE"]:
                continue
            spectra.append(_load_BOSS_HDU(hdulist, hdu, **kwargs))
        return SpectrumList(spectra)


def _load_BOSS_HDU(hdulist: HDUList, hdu: int, **kwargs):
    """
    HDU processor for BOSS spectra redux HDU's

    Parameters
    ----------
    hdulist: HDUList
        HDUList generated from imported file.
    hdu : int
        Specified HDU to load.

    Returns
    -------
    Spectrum1D
        The spectrum contained in the file

    """
    # Fetch flux, e_flux, spectral axis, and units
    # NOTE: flux unit info is not in BOSS spec files.
    flux_unit = Unit("1e-17 erg / (Angstrom cm2 s)")  # NOTE: hardcoded unit
    spectral_axis = Quantity(10**hdulist[hdu].data["LOGLAM"], unit=Angstrom)

    flux = Quantity(hdulist[hdu].data["FLUX"], unit=flux_unit)
    # no e_flux, so we use inverse of variance
    ivar = InverseVariance(hdulist[hdu].data["IVAR"])

    # Fetch mask -- stored in two ways
    try:
        mask = hdulist[hdu].data["MASK"]
    except KeyError:
        mask = hdulist[hdu].data["OR_MASK"]

    # NOTE: specutils considers 0/False as valid values, simlar to numpy convention
    mask = mask != 0

    # Fetch metadata
    # NOTE: specFull file does not include S/N value, but this gets calculated
    #       for mwmVisit/mwmStar files when they are created.
    meta = dict()
    meta["header"] = hdulist[0].header
    meta["name"] = hdulist[hdu].name

    return Spectrum1D(spectral_axis=spectral_axis,
                      flux=flux,
                      uncertainty=ivar,
                      mask=mask,
                      meta=meta)


# MWM LOADERS
@data_loader(
    "SDSS-V mwm",
    identifier=mwm_identify,
    dtype=Spectrum1D,
    priority=20,
    extensions=["fits"],
)
def load_sdss_mwm_1d(file_obj, hdu: Optional[int] = None, **kwargs):
    """
    Load an unspecified spec file as a Spectrum1D.

    Parameters
    ----------
    file_obj : str, file-like, or HDUList
        FITS file name, file object, or HDUList..
    hdu : int
        Specified HDU to load.

    Returns
    -------
    Spectrum1D
        The spectrum contained in the file
    """
    with read_fileobj_or_hdulist(file_obj, memmap=False, **kwargs) as hdulist:
        # Check if file is empty first
        datasums = []
        for i in range(1, len(hdulist)):
            datasums.append(int(hdulist[i].header.get("DATASUM")))
        if (np.array(datasums) == 0).all():
            raise ValueError("Specified file is empty.")

        # TODO: how should we handle this -- multiple things in file, but the user cannot choose.
        if hdu is None:
            for i in range(len(hdulist)):
                if hdulist[i].header.get("DATASUM") != "0":
                    hdu = i
                    break

        return _load_mwmVisit_or_mwmStar_hdu(hdulist, hdu, **kwargs)


@data_loader(
    "SDSS-V mwm multi",
    identifier=mwm_identify,
    dtype=SpectrumList,
    priority=20,
    extensions=["fits"],
)
def load_sdss_mwm_list(file_obj, **kwargs):
    """
    Load an mwmStar spec file as a SpectrumList.

    Parameters
    ----------
    file_obj : str, file-like, or HDUList
        FITS file name, file object, or HDUList.

    Returns
    -------
    SpectrumList
        The spectra contained in the file, where:
            Spectrum1D
                A given spectra of nD flux
            None
                If there are no spectra for that spectrograph/observatory
    """
    spectra = SpectrumList()
    with read_fileobj_or_hdulist(file_obj, memmap=False, **kwargs) as hdulist:
        # Check if file is empty first
        datasums = []
        for hdu in range(1, len(hdulist)):
            datasums.append(int(hdulist[hdu].header.get("DATASUM")))
        if (np.array(datasums) == 0).all():
            raise ValueError("Specified file is empty.")

        # Now load file
        for hdu in range(1, len(hdulist)):
            if hdulist[hdu].header.get("DATASUM") == "0":
                # Skip zero data HDU's
                # TODO: validate if we want this printed warning or not.
                # it might get annoying & fill logs with useless alerts.
                print("WARNING: HDU{} ({}) is empty.".format(
                    hdu, hdulist[hdu].name))
                continue
            spectra.append(_load_mwmVisit_or_mwmStar_hdu(hdulist, hdu))
    return spectra


def _load_mwmVisit_or_mwmStar_hdu(hdulist: HDUList, hdu: int, **kwargs):
    """
    HDU loader subfunction for MWM files

    Parameters
    ----------
    hdulist: HDUList
        HDUList generated from imported file.
    hdu: int
        Specified HDU to load.

    Returns
    -------
    Spectrum1D
        The spectrum with nD flux contained in the HDU.

    """
    if hdulist[hdu].header.get("DATASUM") == "0":
        raise IndexError(
            "Attemped to load an empty HDU specified at HDU{}".format(hdu))

    # Fetch wavelength
    # encoded as WCS for visit, and 'wavelength' for star
    try:
        wavelength = np.array(hdulist[hdu].data["wavelength"])[0]
    except KeyError:
        wavelength = _wcs_log_linear(
            hdulist[hdu].header.get("NPIXELS"),
            hdulist[hdu].header.get("CDELT"),
            hdulist[hdu].header.get("CRVAL"),
        )
    finally:
        if wavelength is None:
            raise ValueError(
                "Couldn't find wavelength data in HDU{}.".format(hdu))
        spectral_axis = Quantity(wavelength, unit=Angstrom)

    # Fetch flux, e_flux
    # NOTE:: flux info is not written into mwm files
    flux_unit = Unit("1e-17 erg / (Angstrom cm2 s)")  # NOTE: hardcoded unit
    flux = Quantity(hdulist[hdu].data["flux"], unit=flux_unit)
    e_flux = InverseVariance(array=hdulist[hdu].data["ivar"])

    # Collect bitmask
    mask = hdulist[hdu].data["pixel_flags"]
    # NOTE: specutils considers 0/False as valid values, simlar to numpy convention
    mask = mask != 0

    # collapse shape if 1D spectra in 2D array
    # NOTE: this fixes a jdaviz handling bug for 2D of shape 1,
    #       it could be that it's expected to be parsed this way.
    if flux.shape[0] == 1:
        flux = flux[0]
        e_flux = e_flux[0]
        mask = mask[0]

    # Create metadata
    meta = dict()
    meta["header"] = hdulist[0].header

    # Add SNR to metadata
    meta["snr"] = np.array(hdulist[hdu].data["snr"])

    # Add identifiers (obj, telescope, mjd, datatype)
    # TODO: need to see what metadata we're interested in for the MWM files.
    meta["telescope"] = hdulist[hdu].data["telescope"]
    meta["instrument"] = hdulist[hdu].header.get("INSTRMNT")
    try:
        meta["obj"] = hdulist[hdu].data["obj"]
    except KeyError:
        pass
    try:
        meta["date"] = hdulist[hdu].data["date_obs"]
        meta["mjd"] = hdulist[hdu].data["mjd"]
        meta["datatype"] = "mwmVisit"
    except KeyError:
        meta["mjd"] = (str(hdulist[hdu].data["min_mjd"][0]) + "->" +
                       str(hdulist[hdu].data["max_mjd"][0]))
        meta["datatype"] = "mwmStar"
    finally:
        meta["name"] = hdulist[hdu].name

    return Spectrum1D(
        spectral_axis=spectral_axis,
        flux=flux,
        uncertainty=e_flux,
        mask=mask,
        meta=meta,
    )
