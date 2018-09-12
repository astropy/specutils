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

    left_index = int(np.ceil(spectrum.wcs.world_to_pixel([subregion[0]])))
    right_index = int(np.floor(spectrum.wcs.world_to_pixel([subregion[1]]))) + 1

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

    """

    extracted_spectrum = []
    for subregion in region._subregions:
        left_index, right_index = to_pixel(subregion, spectrum)
        extracted_spectrum.append(spectrum[left_index:right_index])

    # If there is only one subregion in the region then we will
    # just return a spectrum.
    if len(region) == 1:
        extracted_spectrum = extracted_spectrum[0]

    return extracted_spectrum
