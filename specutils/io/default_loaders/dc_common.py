from collections.abc import Callable
from copy import deepcopy
from enum import Enum

from astropy.nddata import (
    VarianceUncertainty,
    StdDevUncertainty,
    InverseVariance,
)
import astropy.units as u
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_pixel

from ..parsing_utils import read_fileobj_or_hdulist
from ..registers import data_loader
from ... import Spectrum1D, SpectrumList


HEADER_PUPOSE_KEYWORDS = ["EXTNAME", "HDUNAME"]
HEADER_INDEX_PUPOSE_KEYWORDS = ["ROW", "ARRAY"]
FITS_FILE_EXTS = ["fit", "fits", "fts"]
SINGLE_SPLIT_LABEL = "Data Central Single-Split"
MULTILINE_SINGLE_LABEL = "Data Central Multiline-Single"
UNKNOWN_LABEL = "Unable to find a sensible label for spectrum"
# These are order in a best guess of priority, ideally the loader would know
# which label to use.
HEADER_LABEL_KEYWORDS = [
    "OBJECT",
    "OBJNAME",
    "OBS_ID",
    "EXTNAME",
    "HDUNAME",
    "TITLE",
    "ORIGIN",
    "ROOTNAME",
    "FILENAME",
    "AUTHOR",
    "OBSERVER",
    "CREATOR",
    "INSTRUME",
    "PROGRAM",
]


def guess_label_from_header(header):
    """
    Guess the label from `header`, which is assumed to be some mapping with
    FITS-like keys.
    """
    for header_key in HEADER_LABEL_KEYWORDS:
        label = header.get(header_key)
        if label is not None:
            return str(label)
    raise ValueError(UNKNOWN_LABEL)


class Purpose(Enum):
    SKIP = "skip"
    SCIENCE = "science"
    ERROR_STDEV = "error_stdev"
    ERROR_VARIANCE = "error_variance"
    ERROR_INVERSEVARIANCE = "error_inversevariance"
    SKY = "sky"
    COMBINED_SCIENCE = "combined_science"
    COMBINED_ERROR_STDEV = "combined_error_stdev"
    COMBINED_ERROR_VARIANCE = "combined_error_variance"
    COMBINED_ERROR_INVERSEVARIANCE = "combined_error_inversevariance"
    UNREDUCED_SCIENCE = "unreduced_science"
    UNREDUCED_ERROR_STDEV = "unreduced_error_stdev"
    UNREDUCED_ERROR_VARIANCE = "unreduced_error_variance"
    UNREDUCED_ERROR_INVERSEVARIANCE = "unreduced_error_inversevariance"
    NORMALISED_SCIENCE = "normalised_science"
    NORMALISED_ERROR_STDEV = "normalised_error_stdev"
    NORMALISED_ERROR_VARIANCE = "normalised_error_variance"
    NORMALISED_ERROR_INVERSEVARIANCE = "normalised_error_inversevariance"


CREATE_SPECTRA = {
    Purpose.SCIENCE,
    Purpose.SKY,
    Purpose.COMBINED_SCIENCE,
    Purpose.UNREDUCED_SCIENCE,
    Purpose.NORMALISED_SCIENCE,
}
ERROR_PURPOSES = {
    Purpose.ERROR_STDEV,
    Purpose.ERROR_VARIANCE,
    Purpose.ERROR_INVERSEVARIANCE,
    Purpose.COMBINED_ERROR_STDEV,
    Purpose.COMBINED_ERROR_VARIANCE,
    Purpose.COMBINED_ERROR_INVERSEVARIANCE,
    Purpose.UNREDUCED_ERROR_STDEV,
    Purpose.UNREDUCED_ERROR_VARIANCE,
    Purpose.UNREDUCED_ERROR_INVERSEVARIANCE,
    Purpose.NORMALISED_ERROR_STDEV,
    Purpose.NORMALISED_ERROR_VARIANCE,
    Purpose.NORMALISED_ERROR_INVERSEVARIANCE,
}
PURPOSE_SPECTRA_MAP = {
    Purpose.SCIENCE: "reduced",
    Purpose.ERROR_STDEV: "reduced",
    Purpose.ERROR_VARIANCE: "reduced",
    Purpose.ERROR_INVERSEVARIANCE: "reduced",
    Purpose.SKY: "sky",
    Purpose.COMBINED_SCIENCE: "combined",
    Purpose.COMBINED_ERROR_STDEV: "combined",
    Purpose.COMBINED_ERROR_VARIANCE: "combined",
    Purpose.COMBINED_ERROR_INVERSEVARIANCE: "combined",
    Purpose.UNREDUCED_SCIENCE: "unreduced",
    Purpose.UNREDUCED_ERROR_STDEV: "unreduced",
    Purpose.UNREDUCED_ERROR_VARIANCE: "unreduced",
    Purpose.UNREDUCED_ERROR_INVERSEVARIANCE: "unreduced",
    Purpose.NORMALISED_SCIENCE: "normalised",
    Purpose.NORMALISED_ERROR_STDEV: "normalised",
    Purpose.NORMALISED_ERROR_VARIANCE: "normalised",
    Purpose.NORMALISED_ERROR_INVERSEVARIANCE: "normalised",
}
UNCERTAINTY_MAP = {
    Purpose.ERROR_STDEV: StdDevUncertainty,
    Purpose.ERROR_VARIANCE: VarianceUncertainty,
    Purpose.ERROR_INVERSEVARIANCE: InverseVariance,
    Purpose.COMBINED_ERROR_STDEV: StdDevUncertainty,
    Purpose.COMBINED_ERROR_VARIANCE: VarianceUncertainty,
    Purpose.COMBINED_ERROR_INVERSEVARIANCE: InverseVariance,
    Purpose.UNREDUCED_ERROR_STDEV: StdDevUncertainty,
    Purpose.UNREDUCED_ERROR_VARIANCE: VarianceUncertainty,
    Purpose.UNREDUCED_ERROR_INVERSEVARIANCE: InverseVariance,
    Purpose.NORMALISED_ERROR_STDEV: StdDevUncertainty,
    Purpose.NORMALISED_ERROR_VARIANCE: VarianceUncertainty,
    Purpose.NORMALISED_ERROR_INVERSEVARIANCE: InverseVariance,
}
GUESS_TO_PURPOSE = {
    "badpix": Purpose.SKIP,
    "": Purpose.SKIP,
    "sky": Purpose.SKY,
    "stdev": Purpose.ERROR_STDEV,
    "sigma": Purpose.ERROR_STDEV,
    "variance": Purpose.ERROR_VARIANCE,
    "spectrum": Purpose.SCIENCE,
}


def refresh_units(wcs):
    """
    Reparse unit strings to ensure aliases have been applied.

    This is needed to handle legacy spellings that are used by existing surveys.
    """
    new_units = [u.Unit(str(cu)) for cu in wcs.wcs.cunit]
    wcs.wcs.cunit = new_units
    return wcs


def add_labels(spec_list, use_purpose=False):
    not_labeled = 0
    label_set = set()
    for spec in spec_list:
        meta = spec.meta
        purpose = meta.get("purpose")
        if use_purpose:
            tail = " (" + str(purpose) + ")"
        else:
            tail = ""
        try:
            meta["label"] = guess_label_from_header(meta["header"]) + tail
        except ValueError:
            not_labeled += 1
        else:
            label_set.add(meta["label"])

    if len(label_set) + not_labeled < len(spec_list):
        # This implies there are duplicates
        for i, spec in enumerate(spec_list, start=1):
            label = spec.meta.get("label")
            if label is not None:
                spec.meta["label"] = label + " #" + str(i)


def compute_wcs_from_keys_and_values(
    header=None,
    *,
    wavelength_unit_keyword=None,
    wavelength_unit=None,
    pixel_reference_point_keyword=None,
    pixel_reference_point=None,
    pixel_reference_point_value_keyword=None,
    pixel_reference_point_value=None,
    pixel_width_keyword=None,
    pixel_width=None,
):
    if wavelength_unit is None:
        if wavelength_unit_keyword is None:
            raise ValueError(
                "Either wavelength_unit or wavelength_unit_keyword must be "
                "provided"
            )
        wavelength_unit = u.Unit(header[wavelength_unit_keyword])
    if pixel_reference_point is None:
        if pixel_reference_point_keyword is None:
            raise ValueError(
                "Either pixel_reference_point or "
                "pixel_reference_point_keyword must be provided"
            )
        pixel_reference_point = header[pixel_reference_point_keyword]
    if pixel_reference_point_value is None:
        if pixel_reference_point_value_keyword is None:
            raise ValueError(
                "Either pixel_reference_point_value or "
                "pixel_reference_point_value_keyword must be provided"
            )
        pixel_reference_point_value = header[
            pixel_reference_point_value_keyword
        ]
    if pixel_width is None:
        if pixel_width_keyword is None:
            raise ValueError(
                "Either pixel_width or pixel_width_keyword must be provided"
            )
        pixel_width = header[pixel_width_keyword]

    w = WCS(naxis=1)
    w.wcs.crpix[0] = pixel_reference_point
    w.wcs.crval[0] = pixel_reference_point_value
    w.wcs.cdelt[0] = pixel_width
    w.wcs.cunit[0] = wavelength_unit
    return w


def get_flux_units_from_keys_and_values(
    header,
    *,
    flux_unit_keyword="BUNIT",
    flux_unit=None,
    flux_scale_keyword="BSCALE",
    flux_scale=None,
):
    if flux_unit is None and flux_unit_keyword is None:
        raise ValueError(
            "Either flux_unit or flux_unit_keyword must be provided"
        )
    flux_unit_from_header = header.get(flux_unit_keyword)
    if flux_unit is None and flux_unit_from_header is None:
        raise ValueError(
            "No units found for flux, check flux_unit and flux_unit_keyword"
        )
    flux_unit = u.Unit(flux_unit_from_header or flux_unit)

    flux_scale_from_header = header.get(flux_scale_keyword)
    if flux_scale is None and flux_scale_from_header is None:
        flux_scale = 1
    else:
        flux_scale = flux_scale_from_header or flux_scale
    return flux_scale * flux_unit


def add_single_spectra_to_map(
    spectra_map,
    *,
    header,
    data,
    spec_info=None,
    wcs_info=None,
    units_info=None,
    purpose_prefix=None,
    all_standard_units,
    all_keywords,
    valid_wcs,
    index=None,
    drop_wcs_axes=None,
    fallback_header=None,
):
    spec_wcs_info = {}
    spec_units_info = {}
    if wcs_info is not None:
        spec_wcs_info.update(wcs_info)
    if units_info is not None:
        spec_units_info.update(units_info)

    if spec_info is not None:
        spec_wcs_info.update(spec_info.get("wcs", {}))
        spec_units_info.update(spec_info.get("units", {}))
        purpose = spec_info.get("purpose")
    else:
        purpose = None

    try:
        purpose = get_purpose(
            header,
            purpose=purpose,
            purpose_prefix=purpose_prefix,
            all_keywords=all_keywords,
            index=index,
        )
    except ValueError:
        if fallback_header is None:
            raise
        purpose = get_purpose(
            fallback_header,
            purpose=purpose,
            purpose_prefix=purpose_prefix,
            all_keywords=all_keywords,
            index=index,
        )

    if purpose == Purpose.SKIP:
        return None

    if valid_wcs or not spec_wcs_info:
        try:
            wcs = WCS(header)
        except (ValueError, KeyError):
            if fallback_header is None:
                raise
            wcs = WCS(fallback_header)
        if drop_wcs_axes is not None:
            if isinstance(drop_wcs_axes, Callable):
                wcs = drop_wcs_axes(wcs)
            else:
                wcs = wcs.dropaxis(drop_wcs_axes)
    else:
        try:
            wcs = compute_wcs_from_keys_and_values(header, **spec_wcs_info)
        except ValueError:
            if fallback_header is None:
                raise
            wcs = compute_wcs_from_keys_and_values(
                fallback_header, **spec_wcs_info
            )

    wcs = refresh_units(wcs)

    if all_standard_units:
        spec_units_info = {}
    try:
        flux_unit = get_flux_units_from_keys_and_values(
            header, **spec_units_info
        )
    except ValueError:
        if fallback_header is None:
            raise
        flux_unit = get_flux_units_from_keys_and_values(
            fallback_header, **spec_units_info
        )
    flux = data * flux_unit

    meta = {"header": header, "purpose": PURPOSE_SPECTRA_MAP[purpose]}

    if purpose in CREATE_SPECTRA:
        spectrum = Spectrum1D(wcs=wcs, flux=flux, meta=meta)
        spectra_map[PURPOSE_SPECTRA_MAP[purpose]].append(spectrum)
    elif purpose in ERROR_PURPOSES:
        try:
            spectrum = spectra_map[PURPOSE_SPECTRA_MAP[purpose]][-1]
        except IndexError:
            raise ValueError(f"No spectra to associate with {purpose}")
        aligned_flux = pixel_to_pixel(wcs, spectrum.wcs, flux)
        spectrum.uncertainty = UNCERTAINTY_MAP[purpose](aligned_flux)
        spectrum.meta["uncertainty_header"] = header

    return spectrum


def get_purpose(
    header, *, purpose=None, purpose_prefix=None, all_keywords, index=None
):
    def guess_purpose(header):
        for keyword in HEADER_PUPOSE_KEYWORDS:
            guess = header.get(keyword)
            if guess is not None:
                return GUESS_TO_PURPOSE[guess.strip().lower()]
        return None

    def guess_index_purpose(header, index):
        for keyword in HEADER_INDEX_PUPOSE_KEYWORDS:
            guess = header.get(keyword + str(index))
            if guess is not None:
                return GUESS_TO_PURPOSE[guess.strip().lower()]
        return None

    if all_keywords:
        if index is None:
            guessed_purpose = guess_purpose(header)
            if guessed_purpose is not None:
                return guessed_purpose
            if "XTENSION" not in header:
                # we have a primary HDU, assume science
                return Purpose.SCIENCE
            raise ValueError(
                "Cannot identify purpose, cannot use all_keywords"
            )
        guessed_purpose = guess_index_purpose(header, index)
        if guessed_purpose is not None:
            return guessed_purpose
        raise ValueError("Cannot identify purpose, cannot use all_keywords")
    if purpose is not None:
        return Purpose(purpose)
    if purpose_prefix is not None:
        if index is None:
            return Purpose(header.get(purpose_prefix))
        return Purpose(header.get(purpose_prefix + str(index)))
    raise ValueError(
        "Either all_keywords must be True, or one of purpose or "
        "purpose_prefix must not be None."
    )


def no_auto_identify(*args, **kwargs):
    return False


@data_loader(
    label=SINGLE_SPLIT_LABEL, extensions=FITS_FILE_EXTS, dtype=SpectrumList,
    identifier=no_auto_identify,
)
def load_single_split_file(
    filename,
    *,
    hdus,
    wcs,
    units,
    all_standard_units,
    all_keywords,
    valid_wcs,
    label=True,
    drop_wcs_axes=None,
    fallback_header=None,
):
    spectra_map = {
        "sky": [],
        "combined": [],
        "unreduced": [],
        "normalised": [],
        "reduced": [],
    }

    with read_fileobj_or_hdulist(filename) as fits_file:
        if fallback_header is not None:
            if fallback_header is True:
                fallback_header = 0
            fallback_header = fits_file[fallback_header].header

        hdus = deepcopy(hdus)
        if hdus is not None:
            # extract hdu information and validate it
            cycle = hdus.pop("cycle", None)
            cycle_start = hdus.pop("cycle_start", None)
            purpose_prefix = hdus.pop("purpose_prefix", None)
            if len(hdus) != 0 and cycle is None:
                if len(hdus) < len(fits_file):
                    raise ValueError("Not all HDUs have been specified")
                if len(hdus) > len(fits_file):
                    raise ValueError("Too many HDUs have been specified")
            if cycle is not None and cycle_start is None:
                if len(hdus) == 0:
                    raise ValueError(
                        "If HDUs are not specified, cycle_start must be used"
                    )
                cycle_start = len(hdus)
            if cycle is not None:
                cycle_purpose_prefix = cycle.pop("purpose_prefix", None)
                cycle_length = len(cycle)

                # validate cycle
                if (len(fits_file) - cycle_start) % cycle_length != 0:
                    raise ValueError(
                        "Full cycle cannot be read from fits file"
                    )
        else:
            cycle = None
            cycle_start = 0
            purpose_prefix = None
            cycle_purpose_prefix = None
            cycle_length = 0

        for i, fits_hdu in enumerate(fits_file):
            if cycle is not None and i >= cycle_start:
                hdu_info = cycle.get(str((i - cycle_start) % cycle_length))
                hdu_purpose_prefix = cycle_purpose_prefix
            elif hdus is not None:
                hdu_info = hdus.get(str(i))
                hdu_purpose_prefix = purpose_prefix
            else:
                hdu_info = None
                hdu_purpose_prefix = None

            add_single_spectra_to_map(
                spectra_map,
                data=fits_hdu.data,
                header=fits_hdu.header,
                spec_info=hdu_info,
                wcs_info=wcs,
                units_info=units,
                purpose_prefix=hdu_purpose_prefix,
                all_standard_units=all_standard_units,
                all_keywords=all_keywords,
                valid_wcs=valid_wcs,
                drop_wcs_axes=drop_wcs_axes,
                fallback_header=fallback_header,
            )

    if label:
        add_labels(spectra_map["combined"])
        add_labels(spectra_map["reduced"])
        add_labels(spectra_map["normalised"])
        add_labels(spectra_map["unreduced"], use_purpose=True)
        add_labels(spectra_map["sky"], use_purpose=True)

    return SpectrumList(
        spectra_map["combined"] +
        spectra_map["normalised"] +
        spectra_map["reduced"] +
        spectra_map["unreduced"] +
        spectra_map["sky"]
    )


@data_loader(
    label=MULTILINE_SINGLE_LABEL, extensions=FITS_FILE_EXTS,
    dtype=SpectrumList, identifier=no_auto_identify,
)
def load_multiline_single_file(
    filename,
    *,
    hdu,
    wcs,
    units,
    all_standard_units,
    all_keywords,
    valid_wcs,
    label=True,
    drop_wcs_axes=1,
):
    spectra_map = {
        "sky": [],
        "combined": [],
        "unreduced": [],
        "normalised": [],
        "reduced": [],
    }

    with read_fileobj_or_hdulist(filename) as fits_file:
        fits_header = fits_file[0].header
        fits_data = fits_file[0].data
        hdu = deepcopy(hdu)
        if hdu is not None:
            # extract hdu information and validate it
            if hdu.pop("require_transpose", False):
                fits_data = fits_data.T
            purpose_prefix = hdu.pop("purpose_prefix", None)
            num_rows = fits_data.shape[0]
            if len(hdu) != 0:
                if len(hdu) < num_rows:
                    raise ValueError("Not all rows have been specified")
                if len(hdu) > num_rows:
                    raise ValueError("Too many rows have been specified")
        else:
            purpose_prefix = None

        for i, row in enumerate(fits_data, start=1):
            if hdu is not None:
                row_info = hdu.get(str(i))
            else:
                row_info = None

            add_single_spectra_to_map(
                spectra_map,
                header=fits_header,
                data=row,
                index=i,
                spec_info=row_info,
                wcs_info=wcs,
                units_info=units,
                purpose_prefix=purpose_prefix,
                all_standard_units=all_standard_units,
                all_keywords=all_keywords,
                valid_wcs=valid_wcs,
                drop_wcs_axes=drop_wcs_axes,
            )

    if label:
        add_labels(spectra_map["combined"])
        add_labels(spectra_map["reduced"])
        add_labels(spectra_map["normalised"])
        add_labels(spectra_map["unreduced"], use_purpose=True)
        add_labels(spectra_map["sky"], use_purpose=True)

    return SpectrumList(
        spectra_map["combined"] +
        spectra_map["normalised"] +
        spectra_map["reduced"] +
        spectra_map["unreduced"] +
        spectra_map["sky"]
    )
