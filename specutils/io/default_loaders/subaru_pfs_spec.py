"""
Loader for PFS spectrum files.

https://github.com/Subaru-PFS/datamodel/blob/master/datamodel.txt
"""

import os
import re

import numpy as np
from astropy.nddata import StdDevUncertainty
from astropy.units import Unit

from ...spectra import Spectrum1D
from ..parsing_utils import _fits_identify_by_name, read_fileobj_or_hdulist
from ..registers import data_loader

__all__ = ["identify_pfs_spec", "pfs_spec_loader"]

# This RE matches the file name pattern defined in Subaru-PFS' datamodel.txt :
# "pfsObject-%05d-%05d-%s-%016x-%03d-0x%016x.fits" % (catId, tract, patch, objId,
#                                                     nVisit % 1000, pfsVisitHash)
_spec_pattern = re.compile(
    r"pfsObject-(?P<catId>\d{5})-(?P<tract>\d{5})-(?P<patch>.{3})-(?P<objId>\d{16})-(?P<nVisit>\d{3})-(?P<pfsVisitHash>0x\w{16})\.fits"
)


def identify_pfs_spec(origin, *args, **kwargs):
    """
    Check whether given file is FITS and name matches `_spec_pattern`.
    """

    return _fits_identify_by_name(origin, *args, pattern=_spec_pattern)


@data_loader(
    label="Subaru-pfsObject",
    identifier=identify_pfs_spec,
    extensions=["fits"],
    priority=10,
)
def pfs_spec_loader(file_obj, **kwargs):
    """
    Loader for PFS combined spectrum files.

    Parameters
    ----------
    file_obj : str or file-like
        FITS file name or object (provided from name by Astropy I/O Registry).

    Returns
    -------
    data : Spectrum1D
        The spectrum that is represented by the data in this table.
    """

    # This will fail for file-like objects without 'name' property like `bz2.BZ2File`,
    # workarund needed (or better yet, a scheme to parse the `meta` items from the header).
    if isinstance(file_obj, str):
        file_name = file_obj
    else:
        file_name = file_obj.name

    m = _spec_pattern.match(os.path.basename(file_name))

    with read_fileobj_or_hdulist(file_obj, **kwargs) as hdulist:
        header = hdulist[0].header
        meta = {
            "header": header,
            "tract": m["tract"],
            "patch": m["patch"],
            "catId": m["catId"],
            "objId": m["objId"],
            "nVisit": m["nVisit"],
            "pfsVisitHash": m["pfsVisitHash"],
        }

        # PFS datamodel: https://github.com/Subaru-PFS/datamodel/blob/master/datamodel.txt
        #
        # HDU #8 FLUXTABLE   Binary table  [FITS BINARY TABLE] NOBS*NROW
        # Columns for:
        #   wavelength in units of nm (vacuum)          [64-bit FLOAT]
        #   intensity in units of nJy                   [FLOAT]
        #   intensity error in same units as intensity  [FLOAT]
        #   mask                                        [32-bit INT]
        #
        # In the datamodel, the FLUXTABLE is the 8th HDU, but in fact it's in the 10th HDU with EXTNAME of FLUX_TABLE.
        #
        # NOTE: No backward compatibility is guaranteed for the FLUXTABLE format.
        data = hdulist[10].data["flux"]
        unit = Unit("nJy")

        error = hdulist[10].data["error"]
        uncertainty = StdDevUncertainty(np.sqrt(error))

        wave = hdulist[10].data["wavelength"]
        wave_unit = Unit("nm")

        # NOTE: REFLINE mask is the 15th bit in the mask bits, and can be ignored.
        # TODO: Mask condition needs to be refined once more information becomes available.
        mask = ~np.logical_or(
            # mask==0: good data
            hdulist[10].data["mask"] == 0,
            # mask==2**REFLINE (32768) can be ignored
            hdulist[10].data["mask"] == 2 ** (hdulist[10].header["MP_REFLINE"]),
        )

    return Spectrum1D(
        flux=data * unit,
        spectral_axis=wave * wave_unit,
        uncertainty=uncertainty,
        meta=meta,
        mask=mask,
    )
