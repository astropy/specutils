from specutils import SpectrumList
from specutils.io.registers import data_loader

from .dc_common import (
    FITS_FILE_EXTS, SINGLE_SPLIT_LABEL, MULTILINE_SINGLE_LABEL,
)
from ..parsing_utils import read_fileobj_or_hdulist

GAMA_2QZ_CONFIG = {
    "hdus": None,
    "wcs": {
        "pixel_reference_point_keyword": "CRPIX1",
        "pixel_reference_point_value_keyword": "CRVAL1",
        "pixel_width_keyword": "CD1_1",
        "wavelength_unit": "Angstrom",
    },
    "units": {"flux_unit": "count"},
    "all_standard_units": False,
    "all_keywords": True,
    "valid_wcs": False,
}
GAMA_2SLAQ_QSO_CONFIG = {
    "hdus": None,
    "wcs": {
        "pixel_reference_point_keyword": "CRPIX1",
        "pixel_reference_point_value_keyword": "CRVAL1",
        "pixel_width_keyword": "CDELT1",
        "wavelength_unit": "Angstrom",
    },
    "units": {"flux_unit": "count"},
    "all_standard_units": False,
    "all_keywords": True,
    "valid_wcs": False,
}
GAMA_CONFIG = {
    "hdu": {
        "1": {"purpose": "science", "units": {
            "flux_unit": "10^-17 erg/s/cm^2/A"
        }, },
        "2": {"purpose": "error_stdev"},
        "3": {"purpose": "unreduced_science"},
        "4": {"purpose": "unreduced_error_stdev"},
        "5": {"purpose": "sky"},
    },
    "wcs": None,
    "units": {"flux_unit": "count"},
    "all_standard_units": False,
    "all_keywords": False,
    "valid_wcs": True,
}
GAMA_LT_CONFIG = {
    "hdus": {"0": {"purpose": "science"}, },
    "wcs": {
        "pixel_reference_point_keyword": "CRPIX",
        "pixel_reference_point_value_keyword": "CRVAL",
        "pixel_width_keyword": "CDELT",
        "wavelength_unit": "Angstrom",
    },
    "units": {"flux_unit": "count"},
    "all_standard_units": False,
    "all_keywords": False,
    "valid_wcs": False,
}
MGC_CONFIG = {
    "hdu": None,
    "units": {"flux_unit": "count"},
    "wcs": None,
    "all_standard_units": False,
    "all_keywords": True,
    "valid_wcs": True,
}
WIGGLEZ_CONFIG = {
    "hdus": {
        "0": {"purpose": "science"},
        "1": {"purpose": "error_variance"},
        "2": {"purpose": "skip"},
    },
    "wcs": {
        "pixel_reference_point_keyword": "CRPIX1",
        "pixel_reference_point_value_keyword": "CRVAL1",
        "pixel_width_keyword": "CDELT1",
        "wavelength_unit": "Angstrom",
    },
    "units": {"flux_unit": "count"},
    "all_standard_units": False,
    "all_keywords": False,
    "valid_wcs": False,
}
GAMA_2DFGRS_CONFIG = {
    "hdu": None,
    "wcs": {
        "pixel_reference_point_keyword": "CRPIX1",
        "pixel_reference_point_value_keyword": "CRVAL1",
        "pixel_width_keyword": "CDELT1",
        "wavelength_unit": "Angstrom",
    },
    "units": {"flux_unit": "count"},
    "all_standard_units": False,
    "all_keywords": True,
    "valid_wcs": False,
}


def identify_2qz(origin, *args, **kwargs):
    """
    Identify if the current file is a OzDES file
    """
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        header = hdulist[0].header
        if "2QZ" in header.get("SURVEY", "") and "GAMANAME" in header:
            return True
        return False


@data_loader(
    label="GAMA-2QZ", extensions=FITS_FILE_EXTS, dtype=SpectrumList,
    identifier=identify_2qz, priority=10,
)
def twoqz_loader(filename):
    spectra = SpectrumList.read(
        filename, format=SINGLE_SPLIT_LABEL, **GAMA_2QZ_CONFIG
    )
    return spectra


def identify_2slaq_qso(origin, *args, **kwargs):
    """
    Identify if the current file is a OzDES file
    """
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        header = hdulist[0].header
        if "2SLAQ-QSO" in header.get("SURVEY", "") and "GAMANAME" in header:
            return True
        return False


@data_loader(
    label="GAMA-2SLAQ-QSO", extensions=FITS_FILE_EXTS, dtype=SpectrumList,
    identifier=identify_2slaq_qso, priority=10,
)
def twoslaq_qso_loader(filename):
    spectra = SpectrumList.read(
        filename, format=SINGLE_SPLIT_LABEL, **GAMA_2SLAQ_QSO_CONFIG
    )
    return spectra


def identify_gama(origin, *args, **kwargs):
    """
    Identify if the current file is a OzDES file
    """
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        header = hdulist[0].header
        if "GAMA" in header.get("ORIGIN", "").strip() and "GAMANAME" in header:
            return True
        return False


@data_loader(
    label="GAMA", extensions=FITS_FILE_EXTS, dtype=SpectrumList,
    identifier=identify_gama, priority=10,
)
def gama_loader(filename):
    spectra = SpectrumList.read(
        filename, format=MULTILINE_SINGLE_LABEL, **GAMA_CONFIG
    )
    return spectra


def identify_gama_lt(origin, *args, **kwargs):
    """
    Identify if the current file is a OzDES file
    """
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        header = hdulist[0].header
        if "Liverpool JMU" in header.get("ORIGIN", "") and "GAMANAME" in header:
            return True
        return False


@data_loader(
    label="GAMA-LT", extensions=FITS_FILE_EXTS, dtype=SpectrumList,
    identifier=identify_gama_lt, priority=10,
)
def gama_lt_loader(filename):
    spectra = SpectrumList.read(
        filename, format=SINGLE_SPLIT_LABEL, **GAMA_LT_CONFIG
    )
    return spectra


def identify_mgc(origin, *args, **kwargs):
    """
    Identify if the current file is a OzDES file
    """
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        header = hdulist[0].header
        if "MGC" in header.get("SURVEY", "") and "GAMANAME" in header:
            return True
        return False


@data_loader(
    label="GAMA-MGC", extensions=FITS_FILE_EXTS, dtype=SpectrumList,
    identifier=identify_mgc, priority=10,
)
def mgc_loader(filename):
    spectra = SpectrumList.read(
        filename, format=MULTILINE_SINGLE_LABEL, **MGC_CONFIG
    )
    return spectra


def identify_wigglez(origin, *args, **kwargs):
    """
    Identify if the current file is a OzDES file
    """
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        header = hdulist[0].header
        if "WiggleZ" in header.get("SURVEY", "") and "GAMANAME" in header:
            return True
        return False


@data_loader(
    label="GAMA-WiggleZ", extensions=FITS_FILE_EXTS, dtype=SpectrumList,
    identifier=identify_wigglez, priority=10,
)
def wigglez_loader(filename):
    spectra = SpectrumList.read(
        filename, format=SINGLE_SPLIT_LABEL, **WIGGLEZ_CONFIG
    )
    return spectra


def identify_2dfgrs(origin, *args, **kwargs):
    """
    Identify if the current file is a OzDES file
    """
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        header = hdulist[0].header
        if "2dFGRS" in header.get("SURVEY", "") and "GAMANAME" in header:
            return True
        return False


@data_loader(
    label="GAMA-2dFGRS", extensions=FITS_FILE_EXTS, dtype=SpectrumList,
    identifier=identify_2dfgrs, priority=10,
)
def gama_2dfgrs_loader(filename):
    spectra = SpectrumList.read(
        filename, format=MULTILINE_SINGLE_LABEL, **GAMA_2DFGRS_CONFIG
    )
    return spectra
