from math import floor, ceil  # faster than int(np.floor/ceil(float))

import numpy as np
import warnings

from astropy import units as u
from astropy.wcs import WCS
from gwcs import WCS as GWCS
from ..spectra import Spectrum, SpectralRegion
from ..utils.wcs_utils import gwcs_from_array

__all__ = ['extract_region', 'extract_bounding_spectral_region', 'spectral_slab']


def _edge_value_to_pixel(edge_value, spectrum, order, side, axis=None):
    spectral_axis = spectrum.spectral_axis if axis is None else axis
    try:
        edge_value = edge_value.to(spectral_axis.unit, u.spectral())
    except u.UnitConversionError:
        pass
    if order == 'ascending':
        if side == 'right':
            index = np.searchsorted(spectral_axis, edge_value, side='right')
            if np.isclose(spectral_axis[index-1].value, edge_value.value):
                index += 1
        else:
            index = np.searchsorted(spectral_axis, edge_value, side='left')
        return index

    elif order == 'descending':

        if side == 'left':
            opposite_side = 'right'
        else:
            opposite_side = 'left'

        index = len(spectral_axis) - np.searchsorted(spectral_axis[::-1], edge_value, side=opposite_side)
        return index


def _subregion_to_edge_pixels(subregion, spectrum):
    """
    Calculate and return the left and right indices defined
    by the lower and upper bounds and based on the input
    `~specutils.spectra.spectrum.Spectrum`. The left and right indices will
    be returned.

    Parameters
    ----------
    spectrum: `~specutils.spectra.spectrum.Spectrum`
        The spectrum object from which the region will be extracted.

    Returns
    -------
    left_index, right_index: int, int
        Left and right indices defined by the lower and upper bounds.

    """
    spectral_axis = spectrum.spectral_axis

    # Left/lower side of sub region
    if subregion[0].unit.is_equivalent(u.pix):
        # Pixel handling assumes always ascending
        if not spectral_axis.unit.is_equivalent(u.pix):
            left_index = floor(subregion[0].value)
        else:
            # If the lower bound is larger than the largest value, immediately return nothing
            # Assuming ascending, both bounds are "out of bounds"
            if subregion[0] > spectral_axis[-1]:
                return None, None
            else:
                # Get index of closest value
                # See https://stackoverflow.com/a/26026189 if performance becomes an issue
                left_index = np.nanargmin((np.abs(spectral_axis - subregion[0])))
                # Ensure index is inclusive of region bounds
                if (spectral_axis[left_index] > subregion[0]) and (left_index >= 1):
                    left_index -= 1
    else:
        # Convert lower value to the appropriate axis and compute order on that axis
        try:
            axis_to_use = spectral_axis
            left_reg_in_spec_unit = subregion[0].to(axis_to_use.unit, u.spectral())
        except u.UnitConversionError:
            axis_to_use = _get_axis_in_matching_unit(subregion[0].unit, spectrum)
            left_reg_in_spec_unit = subregion[0].to(axis_to_use.unit, u.spectral())

        order_left = "ascending" if axis_to_use[-1] > axis_to_use[0] else "descending"
        left_index = _edge_value_to_pixel(left_reg_in_spec_unit, spectrum, order_left, "left", axis=axis_to_use)

    # Right/upper side of sub region
    if subregion[1].unit.is_equivalent(u.pix):
        # Pixel handling assumes always ascending
        if not spectral_axis.unit.is_equivalent(u.pix):
            right_index = ceil(subregion[1].value)
        else:
            # If upper bound is smaller than smallest value, immediately return nothing
            # Assuming ascending, both bounds are "out of bounds"
            if subregion[1] < spectral_axis[0]:
                return None, None
            else:
                # Get index of closest value
                # See https://stackoverflow.com/a/26026189 if performance becomes an issue
                right_index = np.nanargmin((np.abs(spectral_axis - subregion[1])))
                # Ensure index is inclusive of region bounds
                if (spectral_axis[right_index] < subregion[1]) and (right_index < len(spectral_axis)):
                    right_index += 1
    else:
        # Convert upper value to the appropriate axis and compute order on that axis
        try:
            axis_to_use_r = spectral_axis
            right_reg_in_spec_unit = subregion[1].to(axis_to_use_r.unit, u.spectral())
        except u.UnitConversionError:
            axis_to_use_r = _get_axis_in_matching_unit(subregion[1].unit, spectrum)
            right_reg_in_spec_unit = subregion[1].to(axis_to_use_r.unit, u.spectral())

        # Compute stop-exclusive index so [lower, upper] is closed in world coords
        axis_vals = axis_to_use_r.value
        upper_val = right_reg_in_spec_unit.value

        if axis_to_use_r[-1] > axis_to_use_r[0]:
            # ascending: first index strictly greater than upper
            right_index = int(np.searchsorted(axis_vals, upper_val, side="right"))
        else:
            # descending: search on reversed, then map back
            rvals = axis_vals[::-1]
            right_r = int(np.searchsorted(rvals, upper_val, side="right"))
            right_index = len(axis_vals) - right_r

    # If the spectrum is in wavelength and region is in Hz (for example), these still might be reversed
    if left_index < right_index:
        return left_index, right_index
    else:
        return right_index, left_index


def extract_region(spectrum, region, return_single_spectrum=False, preserve_wcs=False):
    """
    Extract a region from the input `~specutils.Spectrum`
    defined by the lower and upper bounds defined by the ``region``
    instance.  The extracted region will be returned as a new
    `~specutils.Spectrum`.

    Parameters
    ----------
    spectrum: `~specutils.Spectrum`
        The spectrum object from which the region will be extracted.

    region: `~specutils.SpectralRegion`
        The spectral region to extract from the original spectrum.

    return_single_spectrum: `bool`
        If ``region`` has multiple sections, whether to return a single spectrum
        instead of multiple `~specutils.Spectrum` objects.  The returned spectrum
        will be a unique, concatenated, spectrum of all sub-regions.

    preserve_wcs: `bool`
        If True, if the WCS is an astropy.wcs.WCS it will be adjusted and retained in
        the output spectrum(s). If False (default) or the input WCS is a GWCS, the original
        WCS will be dropped and replaced by a lookuptable WCS.

    Returns
    -------
    spectrum: `~specutils.Spectrum` or list of `~specutils.Spectrum`
        Excised spectrum, or list of spectra if the input region contained multiple
        subregions and ``return_single_spectrum`` is `False`.

    Notes
    -----
    The region extracted is a discrete subset of the input spectrum. No interpolation is done
    on the left and right side of the spectrum.

    The region is assumed to be a closed interval (as opposed to Python which is open
    on the upper end).  For example:

        Given:
           A ``spectrum`` with spectral_axis of ``[0.1, 0.2, 0.3, 0.4, 0.5, 0.6]*u.um``.

           A ``region`` defined as ``SpectralRegion(0.2*u.um, 0.5*u.um)``

        And we calculate ``sub_spectrum = extract_region(spectrum, region)``, then the ``sub_spectrum``
        spectral axis will be ``[0.2, 0.3, 0.4, 0.5] * u.um``.

    If the ``region`` does not overlap with the ``spectrum`` then an empty Spectrum object
    will be returned.

    """
    extracted_spectrum = []
    for subregion in region._subregions:
        left_index, right_index = _subregion_to_edge_pixels(subregion, spectrum)

        # If both indices are out of bounds then return an empty spectrum
        if left_index == right_index:
            empty_spectrum = Spectrum(spectral_axis=[]*spectrum.spectral_axis.unit,
                                        flux=[]*spectrum.flux.unit)
            extracted_spectrum.append(empty_spectrum)
        else:
            slices = [slice(None),] * len(spectrum.shape)
            slices[spectrum.spectral_axis_index] = slice(left_index, right_index)
            if len(slices) == 1:
                slices = slices[0]
            else:
                slices = tuple(slices)
            sliced = spectrum[slices]

            # Adjust WCS properly
            if preserve_wcs and isinstance(spectrum.wcs, WCS):
                new_wcs = spectrum.wcs.deepcopy()

                # Set CRPIX = 1.0 (FITS convention: reference pixel is 1-indexed)
                new_wcs.wcs.crpix[0] = 1.0

                # Set CRVAL to match the first spectral axis value in the sliced spectrum
                new_wcs.wcs.crval[0] = sliced.spectral_axis[0].to_value(new_wcs.wcs.cunit[0])

                sliced._wcs = new_wcs

            else:
                if preserve_wcs and isinstance(spectrum.wcs, GWCS):
                    warnings.warn("preserve_wcs does not currently work with GWCS, the result"
                                  " will have a spectral lookup table GWCS.")
                sliced._wcs = gwcs_from_array(sliced._spectral_axis,
                                              sliced.flux.shape,
                                              spectral_axis_index=sliced.spectral_axis_index
                                              )

            extracted_spectrum.append(sliced)

    # If there is only one subregion in the region then we will
    # just return a spectrum.
    if len(region) == 1:
        extracted_spectrum = extracted_spectrum[0]

    # Otherwise, if requested to return a single spectrum, we need to combine
    # the Spectrum objects in extracted_spectrum and return a single object.
    elif return_single_spectrum:
        concat_keys = ['flux', 'uncertainty', 'mask']  # spectral_axis handled manually
        copy_keys = ['velocity_convention', 'rest_value', 'meta']

        # NOTE: WCS is intentionally dropped, which will then fallback on lookup

        def _get_joined_value(sps, key, unique_inds=None):
            if key == 'uncertainty':
                # uncertainty cannot be appended directly as its an object,
                # not an array so instead we'll take a copy of the first entry
                # and overwrite the internal array with an appended array
                uncert = sps[0].uncertainty
                if uncert is None:
                    return None
                uncert._array = np.concatenate([sp.uncertainty._array for sp in sps])
                return uncert[unique_inds] if unique_inds is not None else uncert
            elif key in concat_keys or key == 'spectral_axis':
                if getattr(sps[0], key) is None:
                    return None
                concat_arr = np.concatenate([getattr(sp, key) for sp in sps])
                return concat_arr[unique_inds] if unique_inds is not None else concat_arr
            else:
                # all were extracted from the same input spectrum, so we don't
                # need to make sure the properties match
                return getattr(sps[0], key)

        # we'll need to account for removing overlapped regions in the spectral axis,
        # so we'll concatenate that first and track the unique indices
        spectral_axis = _get_joined_value(extracted_spectrum, 'spectral_axis')
        spectral_axis_unique, unique_inds = np.unique(spectral_axis, return_index=True)
        return Spectrum(spectral_axis=spectral_axis_unique,
                          **{key: _get_joined_value(extracted_spectrum, key, unique_inds)
                             for key in concat_keys+copy_keys})

    return extracted_spectrum


def spectral_slab(spectrum, lower, upper):
    """
    Extract a slab from the input `~specutils.Spectrum`
    defined by the lower and upper bounds defined by the ``region``
    instance.  The extracted region will be returned as a new
    `~specutils.Spectrum`.

    Parameters
    ----------
    spectrum: `~specutils.Spectrum`
        The spectrum object from which the region will be extracted.

    lower, upper: `~astropy.units.Quantity`
        The lower and upper bounds of the region to extract
        from the original spectrum.

    Returns
    -------
    spectrum: `~specutils.Spectrum` or list of `~specutils.Spectrum`
        Excised spectrum, or list of spectra if the input region contained multiple
        subregions.

    Notes
    -----
    This is for now just a proxy for function `extract_region`, to ease the
    transition from spectral-cube.

    """
    region = SpectralRegion(lower, upper)

    return extract_region(spectrum, region)


def extract_bounding_spectral_region(spectrum, region):
    """
    Extract the entire bounding region that encompasses all sub-regions
    contained in a multi-sub-region instance of `~specutils.SpectralRegion`.

    In case only one sub-region exists, this method is equivalent to
    `extract_region`.

    Parameters
    ----------
    spectrum: `~specutils.Spectrum`
        The spectrum object from which the region will be extracted.

    region: `~specutils.SpectralRegion`
        The spectral region to extract from the original spectrum, comprised
        of one or more sub-regions.

    Returns
    -------
    spectrum: `~specutils.Spectrum`
        Excised spectrum from the bounding region defined by the set of
        sub-regions in the input ``region`` instance.

    """
    # If there is only one subregion in the region then we will
    # just return a spectrum.
    if len(region) == 1:
        return extract_region(spectrum, region)

    # Look for indices that bound the entire set of sub-regions.
    min_list = [min(sr) for sr in region._subregions]
    max_list = [max(sr) for sr in region._subregions]

    single_region = SpectralRegion(min(min_list), max(max_list))

    return extract_region(spectrum, single_region)


def _get_axis_in_matching_unit(unit, spectrum):
    """
    Return the appropriate spectral axis (wavelength, frequency, or velocity)
    from the input Spectrum object that matches the given unit.

    Parameters
    ----------
    unit : astropy.units.Unit
        The unit to match (e.g., km/s, Hz, micron).

    spectrum : specutils.Spectrum
        The spectrum from which to select the appropriate axis.

    Returns
    -------
    Quantity
        The corresponding axis: one of spectrum.spectral_axis, spectrum.velocity,
        or spectrum.frequency.

    Raises
    ------
    UnitConversionError
        If the unit is not compatible with any of the known spectral axes.
    """
    if unit.is_equivalent(spectrum.spectral_axis.unit):
        return spectrum.spectral_axis
    elif unit.is_equivalent(u.km / u.s):
        return spectrum.velocity
    elif unit.is_equivalent(u.Hz):
        return spectrum.frequency
    else:
        raise u.UnitConversionError(
            f"Cannot convert subregion unit {unit} to any known spectral axis"
        )
