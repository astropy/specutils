from specutils import SpectrumList
from specutils.io.registers import data_loader

from .dc_common import FITS_FILE_EXTS, SINGLE_SPLIT_LABEL
from ..parsing_utils import read_fileobj_or_hdulist

WIGGLEZ_CONFIG = {
    "hdus": None,
    "wcs": None,
    "units": None,
    "all_standard_units": True,
    "all_keywords": True,
    "valid_wcs": True,
}


def identify_wigglez(origin, *args, **kwargs):
    """
    Identify if the current file is a WiggleZ file
    """
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        if '2018MNRAS.474.4151D' in hdulist[0].header.get("REFCODE"):
            return True
    return False


@data_loader(
    label="WiggleZ", extensions=FITS_FILE_EXTS, dtype=SpectrumList,
    identifier=identify_wigglez, priority=10,
)
def wigglez_loader(fname):
    spectra = SpectrumList.read(
        fname, format=SINGLE_SPLIT_LABEL, **WIGGLEZ_CONFIG
    )
    return spectra
