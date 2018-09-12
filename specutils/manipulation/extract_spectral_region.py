import numpy as np

__all__ = ['extract_region']


def to_pixel(subregion, spectrum):
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

    try:
        left_index = int(np.ceil(spectrum.wcs.world_to_pixel([subregion[0]])))
    except Exception as e:
        left_index = None

    try:
        right_index = int(np.floor(spectrum.wcs.world_to_pixel([subregion[1]]))) + 1
    except Exception as e:
        right_index = None

    return left_index, right_index


def extract_region(spectrum, region):
    """
    Extract a region from the input `~specutils.spectra.spectrum1d.Spectrum1D`
    defined by the lower and upper bounds defined by this SpectralRegion
    instance.  The extracted region will be returned.

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

    """

    extracted_spectrum = []
    for subregion in region._subregions:
        left_index, right_index = to_pixel(subregion, spectrum)

        # If both indices are out of bounds then return None
        if left_index is None and right_index is None:
            extracted_spectrum.append(None)
        else:

            # If only one index is out of bounds then set it to
            # the lower or upper extent
            if left_index is None:
                left_index = 0

            if right_index is None:
                right_index = len(spectrum.spectral_axis)

            extracted_spectrum.append(spectrum[left_index:right_index])

    # If there is only one subregion in the region then we will
    # just return a spectrum.
    if len(region) == 1:
        extracted_spectrum = extracted_spectrum[0]

    return extracted_spectrum
