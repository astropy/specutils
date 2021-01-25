from specutils import SpectrumList
from specutils.io.registers import data_loader

from .dc_common import FITS_FILE_EXTS, SINGLE_SPLIT_LABEL
from ..parsing_utils import read_fileobj_or_hdulist

GALAH_5EXT_CONFIG = {
    "hdus": {
        "0": {
            "purpose": "science",
            "units": {"flux_unit": "count"},
        },
        "1": {
            "purpose": "error_stdev",
            "units": {"flux_unit": "count"},
        },
        "2": {
            "purpose": "unreduced_science",
            "units": {"flux_unit": "count"},
        },
        "3": {
            "purpose": "unreduced_error_stdev",
            "units": {"flux_unit": "count"},
        },
        "4": {
            "purpose": "normalised_science",
            "units": {"flux_unit": ""},
        },
    },
    "wcs": {
        "pixel_reference_point_keyword": "CRPIX1",
        "pixel_reference_point_value_keyword": "CRVAL1",
        "pixel_width_keyword": "CDELT1",
        "wavelength_unit": "Angstrom",
    },
    "units": None,
    "all_standard_units": False,
    "all_keywords": False,
    "valid_wcs": False,
}
GALAH_4EXT_CONFIG = {
    "hdus": {
        "0": {"purpose": "science"},
        "1": {"purpose": "error_stdev"},
        "2": {"purpose": "unreduced_science"},
        "3": {"purpose": "unreduced_error_stdev"},
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


def identify_galah(origin, *args, **kwargs):
    """
    Identify if the current file is a GALAH file
    """
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        if "GALAH" in hdulist[0].header.get("ORIGIN"):
            return True
        return False


@data_loader(
    label="GALAH", extensions=FITS_FILE_EXTS, dtype=SpectrumList,
    identifier=identify_galah, priority=10,
)
def galah_loader(filename):
    with read_fileobj_or_hdulist(filename) as hdulist:
        if len(hdulist) == 5:
            spectra = SpectrumList.read(
                hdulist, format=SINGLE_SPLIT_LABEL, **GALAH_5EXT_CONFIG
            )
            spectra[0].meta["galah_hdu_format"] = 5
        elif len(hdulist) == 4:
            spectra = SpectrumList.read(
                hdulist, format=SINGLE_SPLIT_LABEL, **GALAH_4EXT_CONFIG
            )
            spectra[0].meta["galah_hdu_format"] = 4
        else:
            raise RuntimeError(
                "Unknown GALAH format, has {} extensions".format(len(hdulist))
            )
        return spectra
