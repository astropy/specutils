import sys

from math import floor, ceil  # faster than int(np.floor/ceil(float))

import numpy as np

from astropy import units as u
from ..spectra import Spectrum1D, SpectralRegion

__all__ = ['extract_region', 'extract_bounding_spectral_region', 'spectral_slab']


def _to_edge_pixel(subregion, spectrum):
    """
    Calculate and return the left and right indices defined
    by the lower and upper bounds and based on the input
    `~specutils.spectra.spectrum1d.Spectrum1D`. The left and right indices will
    be returned.

    Parameters
    ----------
    spectrum: `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object from which the region will be extracted.

    Returns
    -------
    left_index, right_index: int, int
        Left and right indices defined by the lower and upper bounds.

    """
    # TODO: spectral regions cannot handle strictly ascending spectral axis
    #  values. Instead, convert to length space if axis given in a desceninding
    #  unit space (e.g. frequency). We should find a more elegant solution.
    spectral_axis = spectrum.spectral_axis

    if spectrum.spectral_axis.unit.physical_type != 'length' and \
            spectrum.spectral_axis.unit.is_equivalent(
                u.AA, equivalencies=u.spectral()):
        spectral_axis = spectrum.spectral_axis.to(
            u.AA, equivalencies=u.spectral())

    #
    # Left/lower side of sub region
    #
    if subregion[0].unit.is_equivalent(u.pix):
        left_index = floor(subregion[0].value)
    else:
        # Convert lower value to spectrum spectral_axis units
        left_reg_in_spec_unit = subregion[0].to(spectral_axis.unit,
                                                u.spectral())

        if left_reg_in_spec_unit < spectral_axis[0]:
            left_index = 0
        elif left_reg_in_spec_unit > spectral_axis[-1]:
            left_index = len(spectrum.spectral_axis)
        else:
            try:
                if hasattr(spectrum.wcs, "spectral"):
                    left_index = spectrum.wcs.spectral.world_to_pixel(left_reg_in_spec_unit)
                else:
                    left_index = spectrum.wcs.world_to_pixel(left_reg_in_spec_unit)
                left_index = int(np.ceil(left_index))
            except Exception as e:
                raise ValueError(
                    "Lower value, {}, could not convert using spectrum's WCS "
                    "{}. Exception: {}".format(
                        left_reg_in_spec_unit, spectrum.wcs, e))

    #
    # Right/upper side of sub region
    #
    if subregion[1].unit.is_equivalent(u.pix):
        right_index = ceil(subregion[1].value)
    else:
        # Convert upper value to spectrum spectral_axis units
        right_reg_in_spec_unit = subregion[1].to(spectral_axis.unit,
                                                 u.spectral())

        if right_reg_in_spec_unit > spectral_axis[-1]:
            right_index = len(spectrum.spectral_axis)
        elif right_reg_in_spec_unit < spectral_axis[0]:
            right_index = 0
        else:
            try:
                if hasattr(spectrum.wcs, "spectral"):
                    right_index = spectrum.wcs.spectral.world_to_pixel(right_reg_in_spec_unit)
                else:
                    right_index = spectrum.wcs.world_to_pixel(right_reg_in_spec_unit)
                right_index = int(np.floor(right_index)) + 1
            except Exception as e:
                raise ValueError(
                    "Upper value, {}, could not convert using spectrum's WCS "
                    "{}. Exception: {}".format(
                        right_reg_in_spec_unit, spectrum.wcs, e))

    return left_index, right_index


def extract_region(spectrum, region, return_single_spectrum=False):
    """
    Extract a region from the input `~specutils.Spectrum1D`
    defined by the lower and upper bounds defined by the ``region``
    instance.  The extracted region will be returned as a new
    `~specutils.Spectrum1D`.

    Parameters
    ----------
    spectrum: `~specutils.Spectrum1D`
        The spectrum object from which the region will be extracted.

    region: `~specutils.SpectralRegion`
        The spectral region to extract from the original spectrum.

    return_single_spectrum: `bool`
        If ``region`` has multiple sections, whether to return a single spectrum
        instead of multiple `~specutils.Spectrum1D` objects.  The returned spectrum
        will be a unique, concatenated, spectrum of all sub-regions.

    Returns
    -------
    spectrum: `~specutils.Spectrum1D` or list of `~specutils.Spectrum1D`
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

    If the ``region`` does not overlap with the ``spectrum`` then an empty Spectrum1D object
    will be returned.

    """
    extracted_spectrum = []
    for subregion in region._subregions:
        left_index, right_index = _to_edge_pixel(subregion, spectrum)

        # If both indices are out of bounds then return an empty spectrum
        if left_index is None and right_index is None:
            empty_spectrum = Spectrum1D(spectral_axis=[]*spectrum.spectral_axis.unit,
                                        flux=[]*spectrum.flux.unit)
            extracted_spectrum.append(empty_spectrum)
        else:

            # If only one index is out of bounds then set it to
            # the lower or upper extent
            if left_index is None:
                left_index = 0

            if right_index is None:
                right_index = len(spectrum.spectral_axis)

            if left_index > right_index:
                left_index, right_index = right_index, left_index

            extracted_spectrum.append(spectrum[..., left_index:right_index])

    # If there is only one subregion in the region then we will
    # just return a spectrum.
    if len(region) == 1:
        extracted_spectrum = extracted_spectrum[0]

    # Otherwise, if requested to return a single spectrum, we need to combine
    # the spectrum1d objects in extracted_spectrum and return a single object.
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
        return Spectrum1D(spectral_axis=spectral_axis_unique,
                          **{key: _get_joined_value(extracted_spectrum, key, unique_inds)
                             for key in concat_keys+copy_keys})

    return extracted_spectrum


def spectral_slab(spectrum, lower, upper):
    """
    Extract a slab from the input `~specutils.Spectrum1D`
    defined by the lower and upper bounds defined by the ``region``
    instance.  The extracted region will be returned as a new
    `~specutils.Spectrum1D`.

    Parameters
    ----------
    spectrum: `~specutils.Spectrum1D`
        The spectrum object from which the region will be extracted.

    lower, upper: `~astropy.units.Quantity`
        The lower and upper bounds of the region to extract
        from the original spectrum.

    Returns
    -------
    spectrum: `~specutils.Spectrum1D` or list of `~specutils.Spectrum1D`
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
    spectrum: `~specutils.Spectrum1D`
        The spectrum object from which the region will be extracted.

    region: `~specutils.SpectralRegion`
        The spectral region to extract from the original spectrum, comprised
        of one or more sub-regions.

    Returns
    -------
    spectrum: `~specutils.Spectrum1D`
        Excised spectrum from the bounding region defined by the set of
        sub-regions in the input ``region`` instance.

    """
    # If there is only one subregion in the region then we will
    # just return a spectrum.
    if len(region) == 1:
        return extract_region(spectrum, region)

    min_left = sys.maxsize
    max_right = -sys.maxsize - 1

    # Look for indices that bound the entire set of sub-regions.
    index_list = [_to_edge_pixel(sr, spectrum) for sr in region._subregions]

    for left_index, right_index in index_list:
        if left_index is not None:
            min_left = min(left_index, min_left)
        if right_index is not None:
            max_right = max(right_index, max_right)

    # If both indices are out of bounds then return an empty spectrum
    if min_left is None and max_right is None:
        empty_spectrum = Spectrum1D(spectral_axis=[]*spectrum.spectral_axis.unit,
                                    flux=[]*spectrum.flux.unit)
        return empty_spectrum
    else:
        # If only one index is out of bounds then set it to
        # the lower or upper extent
        if min_left is None:
            min_left = 0

        if max_right is None:
            max_right = len(spectrum.spectral_axis)

        if min_left > max_right:
            min_left, max_right = max_right, min_left

        return spectrum[..., min_left:max_right]
