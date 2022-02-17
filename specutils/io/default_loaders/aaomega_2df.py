from copy import deepcopy

import numpy as np

from astropy.nddata import VarianceUncertainty
from astropy.table import Table
import astropy.units as u
from specutils import Spectrum1D, SpectrumList
from specutils.io.registers import data_loader

from .dc_common import (
    FITS_FILE_EXTS, compute_wcs_from_keys_and_values, add_labels,
)
from ..parsing_utils import read_fileobj_or_hdulist

AAOMEGA_LOADER = "Data Central AAOmega"
AAOMEGA_2DF_FLUX_UNIT = u.Unit("count")
AAOMEGA_2DF_WCS_SETTINGS = {
    "pixel_reference_point_keyword": "CRPIX1",
    "pixel_reference_point_value_keyword": "CRVAL1",
    "pixel_width_keyword": "CDELT1",
    "wavelength_unit": "Angstrom",
}
AAOMEGA_SCIENCE_INDEX = 0
AAOMEGA_FIBRE_INDEX = 2


def identify_aaomega(origin, *args, **kwargs):
    """
    Identify if the current file is a 2dF-AAOmega file reduced by 2dfdr
    """
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        hdu_names = {hdu.name.strip() for hdu in hdulist}
        if "FIBRES" not in hdu_names:
            return False
        if "REDUCTION_ARGS" not in hdu_names:
            return False
        if hdulist[0].header.get("INSTRUME", "").strip() == "AAOMEGA-2dF":
            return True
        return False


@data_loader(
    label=AAOMEGA_LOADER, extensions=FITS_FILE_EXTS, dtype=SpectrumList,
    identifier=identify_aaomega, priority=10,
)
def load_aaomega_file(filename, *args, **kwargs):
    with read_fileobj_or_hdulist(filename, *args, **kwargs) as fits_file:
        fits_header = fits_file[AAOMEGA_SCIENCE_INDEX].header

        # fits_file is the hdulist
        var_idx = None
        rwss_idx = None
        for idx, extn in enumerate(fits_file):
            if extn.name == "VARIANCE":
                var_idx = idx
            if extn.name == "RWSS":
                rwss_idx = idx
        # science data
        fits_data = fits_file[AAOMEGA_SCIENCE_INDEX].data

        # read in Fibre table data....
        ftable = Table(fits_file[AAOMEGA_FIBRE_INDEX].data)

        # A SpectrumList to hold all the Spectrum1D objects
        sl = SpectrumList()

        # the row var contains the pixel data from the science frame
        for i, row in enumerate(fits_data):
            # Definitely need deepcopy here, otherwise it does *NOT* work!
            fib_header = deepcopy(fits_header)
            # Adjusting some values from primary header so individual fibre
            # spectra have meaningful headers
            fib_header["FLDNAME"] = (
                fits_header["OBJECT"],
                "Name of 2dF .fld file"
            )
            fib_header["FLDRA"] = (
                fits_header["MEANRA"],
                "Right Ascension of 2dF field",
            )
            fib_header["FLDDEC"] = (
                fits_header["MEANDEC"],
                "Declination of 2dF field"
            )
            # Now for the fibre specific information from the Fibre Table
            # (extension 2)
            # Caution: RA and DEC are stored in RADIANS in the FIBRE TABLE!
            fib_header["RA"] = (
                ftable["RA"][i] * 180.0 / np.pi,
                "Right Ascension of fibre from configure .fld file",
            )
            fib_header["DEC"] = (
                ftable["DEC"][i] * 180.0 / np.pi,
                "Declination of fibre from configure .fld file",
            )
            fib_header["OBJECT"] = (
                ftable["NAME"][i],
                "Name of target observed by fibre",
            )
            fib_header["OBJCOM"] = (
                ftable["COMMENT"][i],
                "Comment from configure .fld file for target",
            )
            fib_header["OBJMAG"] = (
                ftable["MAGNITUDE"][i],
                "Magnitude of target observed by fibre",
            )
            fib_header["OBJTYPE"] = (
                ftable["TYPE"][i],
                "Type of target observed by fibre",
            )
            fib_header["OBJPIV"] = (
                ftable["PIVOT"][i],
                "Pivot number used to observe target",
            )
            fib_header["OBJPID"] = (
                ftable["PID"][i],
                "Program ID from configure .fld file",
            )
            fib_header["OBJX"] = (
                ftable["X"][i],
                "X coord of target observed by fibre (microns)",
            )
            fib_header["OBJY"] = (
                ftable["Y"][i],
                "Y coord of target observed by fibre (microns)",
            )
            fib_header["OBJXERR"] = (
                ftable["XERR"][i],
                "X coord error of target observed by fibre (microns)",
            )
            fib_header["OBJYERR"] = (
                ftable["YERR"][i],
                "Y coord error of target observed by fibre (microns)",
            )
            fib_header["OBJTHETA"] = (
                ftable["THETA"][i],
                "Angle of fibre used to observe target",
            )
            fib_header["OBJRETR"] = (
                ftable["RETRACTOR"][i],
                "Retractor number used to observe target",
            )
            # WLEN added around 2005 according to AAOmega obs manual...
            # so not always available
            if "WLEN" in ftable.colnames:
                fib_header["OBJWLEN"] = (
                    ftable["WLEN"][i],
                    "Retractor of target observed by fibre",
                )

            # ftable['TYPE'][i]:
            #   P == program (science)
            #   S == sky
            #   U == unallocated or unused
            #   F == fiducial (guide) fibre
            #   N == broken, dead or no fibre
            meta = {"header": fib_header}

            if ftable["TYPE"][i] == "P":
                meta["purpose"] = "reduced"
            elif ftable["TYPE"][i] == "S":
                meta["purpose"] = "sky"
            else:
                # Don't include other fibres that are not science or sky
                continue

            wcs = compute_wcs_from_keys_and_values(
                fib_header, **AAOMEGA_2DF_WCS_SETTINGS
            )
            flux = row * AAOMEGA_2DF_FLUX_UNIT
            meta["fibre_index"] = i

            # Our science spectrum
            spectrum = Spectrum1D(wcs=wcs, flux=flux, meta=meta)
            # If the VARIANCE spectrum exists, add it as an additional spectrum
            # in the meta dict with key 'variance'
            if var_idx is not None:
                var_data = fits_file[var_idx].data
                var_flux = var_data[i] * AAOMEGA_2DF_FLUX_UNIT ** 2
                spectrum.uncertainty = VarianceUncertainty(var_flux)
            # If the RWSS spectrum exists, add it as an additional spectrum in
            # the meta dict with key 'science_sky'
            # This is an optional extension produced by 2dfdr on request: all
            # spectra without the average/median sky subtraction
            # Useful in case users want to do their own sky subtraction.
            if rwss_idx is not None:
                rwss_data = fits_file[rwss_idx].data
                rwss_flux = rwss_data[i] * AAOMEGA_2DF_FLUX_UNIT
                rwss_meta = {
                    "header": fib_header,
                    "purpose": "science_sky"
                }
                spectrum.meta["science_sky"] = Spectrum1D(
                    wcs=wcs, flux=rwss_flux, meta=rwss_meta
                )

            # Add our spectrum to the list.
            # The additional spectra are accessed using
            # spectrum.meta['variance'] and spectrum.meta['science_sky']
            sl.append(spectrum)

    add_labels(sl)

    return sl
