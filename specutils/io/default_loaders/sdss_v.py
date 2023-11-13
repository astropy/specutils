"""Register reader functions for various spectral formats."""
from collections import OrderedDict

import numpy as np
from astropy.units import Unit, Quantity, Angstrom
from astropy.nddata import StdDevUncertainty, InverseVariance

from ...spectra import Spectrum1D, SpectrumList
from ..registers import data_loader
from ..parsing_utils import read_fileobj_or_hdulist

__all__ = ["load_sdss_mwmStar_1d", "load_sdss_mwmVisit_1d"]
"""
From the specutils documentation:
- For spectra that have different shapes, use SpectrumList.
- For spectra that have the same shape but different spectral axes, see SpectrumCollection. 
- For a spectrum or spectra that all share the same spectral axis, use Spectrum1D. 
"""


## IDENTIFIER
def identify_filetype(filetype):
    """Generic identifier function creator for a given filetype based on filepath."""
    return lambda f, o, *args, **kwargs: o.split("/")[-1].startswith(filetype)


## HELPERS
def _wcs_log_linear(naxis, cdelt, crval):
    """
    Convert WCS from log to linear.
    """
    return 10**(np.arange(naxis) * cdelt + crval)


def _fetch_metadata(image):
    """
    Fetch the relevant common metadata.

    Global function used in all SDSS-V metadata handling.

    Parameters
    ----------
    image: HDUList
        The image imported from fits.open().

    Returns
    -------
    meta: OrderedDict
        Dictionary of relevant metadata from primary HDU.

    """
    # Orderly with an OrderedDict
    common_meta = OrderedDict([])

    # Loop over all keys to get relevant metadata
    for key in image[0].header.keys():
        # exclude these keys
        if key.startswith(("TTYPE", "TFORM", "TDIM")) or key in (
                "",
                "COMMENT",
                "CHECKSUM",
                "DATASUM",
                "NAXIS",
                "NAXIS1",
                "NAXIS2",
                "XTENSION",
                "BITPIX",
                "PCOUNT",
                "GCOUNT",
                "TFIELDS",
        ):
            continue

        common_meta[key.lower()] = image[0].header[key]  # add key to dict

    return common_meta


def _fetch_flux_unit(hdu):
    """
    Fetch the flux unit by accessing the BUNIT key of the given HDU.

    Parameters
    ----------
    hdu : flux or e_flux hdu

    Returns
    -------
    flux_unit : astropy Units object of flux unit

    """
    flux_unit = hdu.header["BUNIT"].replace("(",
                                            "").replace(")",
                                                        "").split(" ")[-2:]
    flux_unit = "".join(flux_unit)
    if "Ang" in flux_unit and "strom" not in flux_unit:
        flux_unit = flux_unit.replace("Ang", "Angstrom")

    # TODO: make it so the string is nice, so the stupid warning is supressed
    return Unit(flux_unit)


## APOGEE files
# TODO: add LSF and wavelength coefficients to metadata?
@data_loader(
    "SDSS-V apStar",
    identifier=identify_filetype(("apStar", "asStar")),
    dtype=Spectrum1D,
    priority=10,
    extensions=["fits"],
)
def load_sdss_apStar_1D(path, data_slice=None, **kwargs):
    """
    Load an apStar file as a Spectrum1D.

    Parameters
    ----------
    path: str
        The path to the FITS file

    Returns
    -------
    Spectrum1D
        The spectrum contained in the file

    """
    # intialize slicer
    slicer = slice(*data_slice) if data_slice is not None else slice(None)

    with read_fileobj_or_hdulist(path) as image:
        # Obtain spectral axis from header (HDU[1])
        wavelength = _wcs_log_linear(
            image[1].header["NAXIS1"],
            image[1].header["CDELT1"],
            image[1].header["CRVAL1"],
        )
        spectral_axis = Quantity(wavelength, unit=Angstrom)

        # Obtain flux and e_flux from data (HDU[1] & HDU[2])
        flux_unit = _fetch_flux_unit(image[1])
        flux = Quantity(image[1].data[slicer], unit=flux_unit)
        e_flux = StdDevUncertainty(image[2].data[slicer])

        # Obtain stacked SNR from primary HDU based on no. of visits
        snr = [image[0].header["SNR"]]
        n_visits = image[0].header["NVISITS"]
        if n_visits > 1:
            snr.append(
                snr[0])  # duplicate S/N value for second stacking method
            snr.extend([
                image[0].header[f"SNRVIS{i}"] for i in range(1, 1 + n_visits)
            ])

        meta = _fetch_metadata(image)

        meta["SNR"] = np.array(snr)[slicer]
        meta["BITMASK"] = image[3].data[slicer]

    return Spectrum1D(spectral_axis=spectral_axis,
                      flux=flux,
                      uncertainty=e_flux,
                      meta=meta)


# @data_loader(
#    "SDSS-V apStar multi",
#    identifier=identify_filetype(("apStar", "asStar")),
#    dtype=SpectrumList,
#    priority=10,
#    extensions=["fits"],
# )
# def load_sdss_apStar_list(path, **kwargs):
#    """
#    Load an apStar file as a SpectrumList.
#
#    Parameters
#    ----------
#    path: str
#        The path to the FITS file
#
#    Returns
#    -------
#    SpectrumList
#        The spectra contained in the file
#    """
#    # TODO: this case isn't clear to me -- apStar CAN contain resampled
#    # visit spectra, but where in the HDUList?
#    raise NotImplementedError
#    return SpectrumList([load_sdss_apStar_1D(path, **kwargs)])


@data_loader(
    "SDSS-V apVisit",
    identifier=identify_filetype("apVisit"),
    dtype=Spectrum1D,
    priority=10,
    extensions=["fits"],
)
def load_sdss_apVisit_1D(path, **kwargs):
    """
    Load an apVisit file as a Spectrum1D.

    Parameters
    ----------
    path: str
        The path to the FITS file

    Returns
    -------
    Spectrum1D
        The spectrum contained in the file
    """
    flux_unit = Unit("1e-17 erg / (Angstrom cm2 s)")

    with read_fileobj_or_hdulist(path) as image:
        # Fetch and flatten spectrum parameters
        spectral_axis = Quantity(image[4].data.flatten(), unit=Angstrom)
        flux_unit = _fetch_flux_unit(image[1])
        flux = Quantity(image[1].data.flatten(), unit=flux_unit)
        e_flux = StdDevUncertainty(image[2].data.flatten())

        # Get metadata and attach bitmask and MJD
        meta = _fetch_metadata(image)
        meta["mjd"] = image[0].header["MJD5"]
        meta["date"] = image[0].header["DATE-OBS"]
        meta["bitmask"] = image[3].data.flatten()

    return Spectrum1D(spectral_axis=spectral_axis,
                      flux=flux,
                      uncertainty=e_flux,
                      meta=meta)


@data_loader(
    "SDSS-V apVisit multi",
    identifier=identify_filetype("apVisit"),
    dtype=SpectrumList,
    priority=10,
    extensions=["fits"],
)
def load_sdss_apVisit_list(path, **kwargs):
    """
    Load an apVisit file as a SpectrumList.

    Parameters
    ----------
    path: str
        The path to the FITS file

    Returns
    -------
    SpectrumList
        The spectra from each chip contained in the file
    """
    spectra = SpectrumList()
    with read_fileobj_or_hdulist(path) as image:
        # Get metadata
        flux_unit = _fetch_flux_unit(image[1])
        common_meta = _fetch_metadata(image)

        for chip in range(image[1].data.shape[0]):
            # Fetch spectral axis and flux, and E_flux
            spectral_axis = Quantity(image[4].data[chip], unit=Angstrom)
            flux = Quantity(image[1].data[chip], unit=flux_unit)
            e_flux = StdDevUncertainty(image[2].data[chip])

            # Copy metadata for each, adding chip bitmask and MJD to each
            meta = common_meta.copy()
            meta["bitmask"] = image[3].data[chip]
            meta["mjd"] = image[0].header["MJD5"]
            meta["date"] = image[0].header["DATE-OBS"]

            spectra.append(
                Spectrum1D(
                    spectral_axis=spectral_axis,
                    flux=flux,
                    uncertainty=e_flux,
                    meta=meta,
                ))

    return spectra


## BOSS REDUX products (specLite, specFull, custom coadd files, etc)


@data_loader(
    "SDSS-V specFull",
    identifier=identify_filetype("specFull"),
    dtype=Spectrum1D,
    priority=10,
    extensions=["fits"],
)
def load_sdss_specFull_1D(path, **kwargs):
    """
    Load the coadd spectrum in a specFull file as a Spectrum1D.

    Parameters
    ----------
    path: str
        The path to the FITS file

    Returns
    -------
    Spectrum1D
        The spectrum contained in the file
    """
    return _load_BOSS_spec_1D(path, **kwargs)


@data_loader(
    "SDSS-V specLite",
    identifier=identify_filetype("specLite"),
    dtype=Spectrum1D,
    priority=10,
    extensions=["fits"],
)
def load_sdss_specLite_1D(path, **kwargs):
    """
    Load the coadd spectrum in a specLite file as a Spectrum1D.

    Parameters
    ----------
    path: str
        The path to the FITS file

    Returns
    -------
    Spectrum1D
        The spectrum contained in the file
    """
    return _load_BOSS_spec_1D(path, **kwargs)


@data_loader(
    "SDSS-V specFull multi",
    identifier=identify_filetype("specFull"),
    dtype=SpectrumList,
    priority=10,
    extensions=["fits"],
)
def load_sdss_specFull_list(path, **kwargs):
    """
    Load an specFull file as a SpectrumList.

    Parameters
    ----------
    path: str
        The path to the FITS file

    Returns
    -------
    SpectrumList
        The spectra contained in the file
    """
    return _load_BOSS_spec_list(path, **kwargs)


@data_loader(
    "SDSS-V specLite multi",
    identifier=identify_filetype("specLite"),
    dtype=SpectrumList,
    priority=10,
    extensions=["fits"],
)
def load_sdss_specLite_list(path, **kwargs):
    """
    Load an specLite file as a SpectrumList.

    Parameters
    ----------
    path: str
        The path to the FITS file

    Returns
    -------
    SpectrumList
        The spectra contained in the file
    """
    return _load_BOSS_spec_list(path, **kwargs)


@data_loader(
    "SDSS-V specGeneric",
    identifier=identify_filetype("spec"),
    dtype=SpectrumList,
    priority=50,
    extensions=["fits"],
)
def load_sdss_specOther_1D(path, **kwargs):
    """
    Load an unspecified spec file as a Spectrum1D.

    Parameters
    ----------
    path: str
        The path to the FITS file

    Returns
    -------
    Spectrum1D
        The spectrum contained in the file
    """
    return _load_BOSS_spec_1D(path, **kwargs)


@data_loader(
    "SDSS-V specGeneric multi",
    identifier=identify_filetype("spec"),
    dtype=SpectrumList,
    priority=50,
    extensions=["fits"],
)
def load_sdss_specOther_list(path, **kwargs):
    """
    Load an unspecified spec file as a SpectrumList.

    Parameters
    ----------
    path: str
        The path to the FITS file

    Returns
    -------
    SpectrumList
        The spectra contained in the file
    """
    return _load_BOSS_spec_list(path, **kwargs)


def _load_BOSS_spec_1D(path, **kwargs):
    """
    Load a given BOSS spec coadd file as a Spectrum1D object.
    """
    with read_fileobj_or_hdulist(path) as image:
        # directly load the coadd at HDU1
        return _load_BOSS_HDU(image, 1, **kwargs)


def _load_BOSS_spec_list(path, **kwargs):
    """
    Load a given BOSS spec file as a SpectrumList object.

    """
    with read_fileobj_or_hdulist(path) as image:
        spectra = list()
        for hdu in range(1, len(image)):
            if image[hdu].name in ["SPALL", "ZALL", "ZLINE"]:
                continue
            spectra.append(_load_BOSS_HDU(image, hdu, **kwargs))
        return SpectrumList(spectra)


def _load_BOSS_HDU(image, hdu: int, **kwargs):
    """
    Sub-loader for BOSS spectra redux HDU's
    """
    # Fetch flux, e_flux, spectral axis, and units
    # TODO: check with data team about loc of BUNIT in specFull/Lite files
    flux_unit = Unit("1e-17 erg / (Angstrom cm2 s)")
    spectral_axis = Quantity(10**image[hdu].data["LOGLAM"], unit=Angstrom)

    flux = Quantity(image[hdu].data["FLUX"], unit=flux_unit)
    # no e_flux, so we use inverse of variance
    ivar = InverseVariance(image[hdu].data["IVAR"])

    # Fetch metadata
    # NOTE: specFull file does not include S/N value, but this gets calculated
    #       for mwmVisit/mwmStar files when they are created
    meta = _fetch_metadata(image)
    meta["name"] = image[hdu].name

    return Spectrum1D(spectral_axis=spectral_axis,
                      flux=flux,
                      uncertainty=ivar,
                      meta=meta)


## MWM LOADERS


@data_loader(
    "SDSS-V mwmVisit",
    identifier=identify_filetype("mwmVisit"),
    dtype=Spectrum1D,
    priority=1,
    extensions=["fits"],
)
def load_sdss_mwmVisit_1d(path, hdu: int, **kwargs):
    """
    Load an unspecified spec file as a Spectrum1D.

    Parameters
    ----------
    path: str
        The path to the FITS file

    Returns
    -------
    Spectrum1D
        The spectrum contained in the file
    """
    with read_fileobj_or_hdulist(path) as image:
        return _load_mwmVisit_or_mwmStar_hdu(image, hdu, **kwargs)


@data_loader(
    "SDSS-V mwmStar",
    identifier=identify_filetype("mwmStar"),
    dtype=Spectrum1D,
    priority=20,
    extensions=["fits"],
)
def load_sdss_mwmStar_1d(path, hdu: int, **kwargs):
    """
    Load an unspecified spec file as a Spectrum1D.

    Parameters
    ----------
    path: str
        The path to the FITS file

    Returns
    -------
    Spectrum1D
        The spectrum contained in the file
    """
    with read_fileobj_or_hdulist(path) as image:
        return _load_mwmVisit_or_mwmStar_hdu(image, hdu, **kwargs)


@data_loader(
    "SDSS-V mwmVisit multi",
    identifier=identify_filetype("mwmVisit"),
    dtype=SpectrumList,
    priority=1,
    extensions=["fits"],
)
def load_sdss_mwmVisit_list(path, **kwargs):
    """
    Load an mwmVisit spec file as a SpectrumList.

    Parameters
    ----------
    path: str
        The path to the FITS file

    Returns
    -------
    SpectrumList
        The spectra contained in the file
    """
    return _load_mwmVisit_or_mwmStar_list(path, **kwargs)


@data_loader(
    "SDSS-V mwmStar multi",
    identifier=identify_filetype("mwmStar"),
    dtype=SpectrumList,
    priority=20,
    extensions=["fits"],
)
def load_sdss_mwmStar_list(path, **kwargs):
    """
    Load an mwmStar spec file as a SpectrumList.

    Parameters
    ----------
    path: str
        The path to the FITS file

    Returns
    -------
    SpectrumList
        The spectra contained in the file
    """
    return _load_mwmVisit_or_mwmStar_list(path, **kwargs)


def _load_mwmVisit_or_mwmStar_list(path, **kwargs):
    """
    List loader subfunction for MWM files.
    """
    spectra = SpectrumList()
    with read_fileobj_or_hdulist(path) as image:
        for hdu in range(1, len(image)):
            if image[hdu].header["DATASUM"] == "0":
                # TODO: What should we return? (empty file case)
                spectra.append(None)
                continue
            spectra.append(_load_mwmVisit_or_mwmStar_hdu(image, hdu))
    return spectra


def _load_mwmVisit_or_mwmStar_hdu(image, hdu: int, **kwargs):
    """
    HDU loader subfunction for MWM files
    """
    if image[hdu].header["DATASUM"] == "0":
        raise ValueError(
            "Attemped to load an empty HDU specified as HDU{}".format(hdu))

    # Fetch wavelength
    # encoded as WCS for visit, and 'wavelength' for star
    try:
        wavelength = np.array(image[hdu].data["wavelength"])[0]
    except KeyError:
        wavelength = _wcs_log_linear(
            image[hdu].header["NPIXELS"],
            image[hdu].header["CDELT"],
            image[hdu].header["CRVAL"],
        )
    finally:
        if wavelength is None:
            raise ValueError(
                "Couldn't find wavelength data in HDU{}.".format(hdu))
        spectral_axis = Quantity(wavelength, unit=Angstrom)

    # Fetch flux, e_flux
    # TODO: check with data team about key for BUNIT in mwm files
    flux_unit = Unit("1e-17 erg / (Angstrom cm2 s)")
    flux = Quantity(image[hdu].data["flux"], unit=flux_unit)
    e_flux = InverseVariance(array=image[hdu].data["ivar"])

    # Access metadata
    meta = _fetch_metadata(image)

    # Add SNR to metadata
    meta["snr"] = np.array(image[hdu].data["snr"])

    # Add identifiers (obj, telescope, mjd, datatype)
    # TODO: need to see what we actually want an identifier as
    meta["telescope"] = image[hdu].data["telescope"]
    meta["obj"] = image[hdu].data["obj"]
    try:
        meta["date"] = image[hdu].data["date_obs"]
        print(meta["date"])
        meta["mjd"] = image[hdu].data["mjd"]
        print(meta["mjd"])
        meta["datatype"] = "mwmVisit"
    except KeyError:
        meta["mjd"] = (str(image[hdu].data["min_mjd"][0]) + "->" +
                       str(image[hdu].data["max_mjd"][0]))
        meta["datatype"] = "mwmStar"
    finally:
        meta["name"] = image[hdu].name

    return Spectrum1D(spectral_axis=spectral_axis,
                      flux=flux,
                      uncertainty=e_flux,
                      meta=meta)
