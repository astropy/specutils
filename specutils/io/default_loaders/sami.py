from astropy.nddata import VarianceUncertainty
from astropy.table import QTable
from specutils import SpectrumList
from specutils.io.default_loaders.dc_common import (
    FITS_FILE_EXTS, add_single_spectra_to_map,
)
from specutils.io.parsing_utils import read_fileobj_or_hdulist
from specutils.io.registers import data_loader

# There appears to be nothing which says "this is a SAMI 1D spectra", so guess
# it based on the headers that should be there
SAMI_1D_SPECTRA_HEADER_KEYWORDS = [
    "BUNIT", "CATADEC", "CATARA", "CDELT1", "CRPIX1", "CRVAL1", "CTYPE1",
    "CUNIT1", "DROPFACT", "ELLIP", "GRATID", "IFUPROBE", "KPC_SIZE", "NAME",
    "N_SPAX", "POS_ANG", "PSFALPHA", "PSFBETA", "PSFFWHM", "RADESYS", "RADIUS",
    "RO_GAIN", "RO_NOISE", "STDNAME", "WCSAXES", "Z_TONRY",
]


def identify_sami_cube(origin, *args, **kwargs):
    """
    Identify if the current file is a SAMI cube file
    """
    # TODO check this
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        header = hdulist[0].header
        data = hdulist[0].data
        if "SAMI" in header.get("INSTRUME", "") and len(data.shape) == 3:
            return True
        return False


def identify_sami_1d_spec(origin, *args, **kwargs):
    """
    Identify if the current file is a SAMI 1d spectra file of some kind
    """
    # TODO check this
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        header = hdulist[0].header
        for key in SAMI_1D_SPECTRA_HEADER_KEYWORDS:
            if key not in header:
                return False
        return True


@data_loader(
    label="SAMI-cube", extensions=FITS_FILE_EXTS, dtype=SpectrumList,
    identifier=identify_sami_cube, priority=10,
)
def sami_cube_loader(filename):
    spectra_map = {
        "sky": [],
        "combined": [],
        "unreduced": [],
        "reduced": [],
    }
    primary_header = None

    with read_fileobj_or_hdulist(filename) as hdulist:
        for i, hdu in enumerate(hdulist):
            if i == 0:
                # This is the primary extension, and the one with the
                # science data. The header is fairly complete.
                primary_header = hdu.header
                spec = add_single_spectra_to_map(
                    spectra_map,
                    header=primary_header,
                    data=hdu.data,
                    index=None,
                    all_standard_units=True,
                    all_keywords=True,
                    valid_wcs=True,
                )

            elif "VARIANCE" == hdu.header.get("EXTNAME"):
                # This is the variance extension, and is missing wcs and
                # units.
                uncertainty = VarianceUncertainty(
                    hdu.data, unit=spec.flux.unit ** 2
                )
                spec.uncertainty = uncertainty

            elif "WEIGHT" == hdu.header.get("EXTNAME"):
                # This is the weight extension, and is missing wcs. The
                # units are effectively "normalised" (from 0-1 it seems).
                spec.meta["sami_cube_weight_map"] = hdu.data

            elif "COVAR" == hdu.header.get("EXTNAME"):
                # This is the spatial covariance extension. It's not clear
                # as to how best to expose this, so skipping for now.
                pass

            elif "QC" == hdu.header.get("EXTNAME"):
                # This is the QC extension, and is a binary table. This we
                # add to the metadata.
                spec.meta["sami_QC_table"] = QTable.read(hdu)

            elif "DUST" == hdu.header.get("EXTNAME"):
                # This is the dust extension, and is missing wcs and
                # units. This should likely be represented as an array plus
                # the metadata in the header.
                spec.meta["sami_dust_vector_weights"] = hdu.data

            elif "BIN_MASK" == hdu.header.get("EXTNAME"):
                # This is the bin mask extension, where the value of each
                # pixel indicates the bin to which it belongs. The bin mask
                # is used to construct the binned fluxes and variances in
                # the above two extensions from the default cubes.
                # This is not the same as the aperture spectra mask with the
                # same HDU name.
                spec.meta["sami_bin_mask"] = hdu.data

            else:
                raise NotImplementedError(
                    "Extension is not handled: index {}; name {}".format(
                        i, hdu.header.get("EXTNAME")
                    )
                )

    spectra = SpectrumList(
        spectra_map["combined"] +
        spectra_map["reduced"] +
        spectra_map["unreduced"] +
        spectra_map["sky"]
    )
    return spectra


@data_loader(
    label="SAMI-1d-spec", extensions=FITS_FILE_EXTS, dtype=SpectrumList,
    identifier=identify_sami_1d_spec, priority=10,
)
def sami_1d_spec_loader(filename):
    spectra_map = {
        "sky": [],
        "combined": [],
        "unreduced": [],
        "reduced": [],
    }
    primary_header = None

    with read_fileobj_or_hdulist(filename) as hdulist:
        for i, hdu in enumerate(hdulist):
            if i == 0:
                # This is the primary extension, and the one with the
                # science data. The header is fairly complete.
                primary_header = hdu.header
                spec = add_single_spectra_to_map(
                    spectra_map,
                    header=primary_header,
                    data=hdu.data,
                    index=None,
                    all_standard_units=True,
                    all_keywords=True,
                    valid_wcs=True,
                )

            elif "VARIANCE" == hdu.header.get("EXTNAME"):
                # This is the variance extension, and is missing wcs and
                # units.
                uncertainty = VarianceUncertainty(
                    hdu.data, unit=spec.flux.unit ** 2
                )
                spec.uncertainty = uncertainty

            elif "BIN_MASK" == hdu.header.get("EXTNAME"):
                # Contains the bin mask used to construct the aperture
                # spectra. A 1 indicates a spaxel was included in the
                # aperture, a 0 indicates a spaxel was not included.
                spec.meta["sami_aperture_spectra_mask"] = hdu.data

            else:
                raise NotImplementedError(
                    "Extension is not handled: index {}; name {}".format(
                        i, hdu.header.get("EXTNAME")
                    )
                )

    spectra = SpectrumList(
        spectra_map["combined"] +
        spectra_map["reduced"] +
        spectra_map["unreduced"] +
        spectra_map["sky"]
    )
    return spectra
