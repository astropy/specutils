from math import floor, ceil  # faster than int(np.floor/ceil(float))

import numpy as np

from astropy import units as u
from .. import Spectrum1D

__all__ = ['extract_region']


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
            left_index = len(spectrum.spectral_axis)-1
        else:
            try:
                left_index = int(np.ceil(spectrum.wcs.world_to_pixel(
                    left_reg_in_spec_unit)))
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
            right_index = len(spectrum.spectral_axis)-1
        elif right_reg_in_spec_unit < spectral_axis[0]:
            right_index = 0
        else:
            try:
                right_index = int(np.floor(spectrum.wcs.world_to_pixel(
                    right_reg_in_spec_unit))) + 1
            except Exception as e:
                raise ValueError(
                    "Upper value, {}, could not convert using spectrum's WCS "
                    "{}. Exception: {}".format(
                        right_reg_in_spec_unit, spectrum.wcs, e))

    return left_index, right_index


def extract_region(spectrum, region):
    """
    Extract a region from the input `~specutils.Spectrum1D`
    defined by the lower and upper bounds defined by the ``region``
    instance.  The extracted region will be returned as a new
    `~specutils.Spectrum1D`.

    Parameters
    ----------
    spectrum: `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object from which the region will be extracted.

    Returns
    -------
    spectrum: `~specutils.spectra.spectrum1d.Spectrum1D`
        Excised spectrum.

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

        # If both indices are out of bounds then return None
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

            extracted_spectrum.append(spectrum[left_index:right_index])

    # If there is only one subregion in the region then we will
    # just return a spectrum.
    if len(region) == 1:
        extracted_spectrum = extracted_spectrum[0]

    return extracted_spectrum
