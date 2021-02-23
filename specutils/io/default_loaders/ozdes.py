from specutils import SpectrumList
from specutils.io.registers import data_loader

from .dc_common import FITS_FILE_EXTS, SINGLE_SPLIT_LABEL
from ..parsing_utils import read_fileobj_or_hdulist

OZDES_CONFIG = {
    "hdus": {
        "0": {"purpose": "combined_science"},
        "1": {"purpose": "combined_error_variance"},
        "2": {"purpose": "skip"},
        "cycle": {
            "0": {"purpose": "science"},
            "1": {"purpose": "error_variance"},
            "2": {"purpose": "skip"},
        },
    },
    "units": None,
    "wcs": None,
    "all_standard_units": True,
    "all_keywords": False,
    "valid_wcs": True,
}


def identify_ozdes(origin, *args, **kwargs):
    """
    Identify if the current file is a OzDES file
    """
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        if "ozdes" in hdulist[0].header.get("REFERENC"):
            return True
        return False


@data_loader(
    label="OzDES", extensions=FITS_FILE_EXTS, dtype=SpectrumList,
    identifier=identify_ozdes, priority=10
)
def ozdes_loader(filename):
    spectra = SpectrumList.read(
        filename, format=SINGLE_SPLIT_LABEL, **OZDES_CONFIG
    )
    return spectra
