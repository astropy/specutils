"""
A module for analysis tools focused on determining the width of
spectral features.
"""

import numpy as np
from astropy.stats.funcs import gaussian_sigma_to_fwhm
from astropy.modeling.models import Gaussian1D
from ..manipulation import extract_region
from . import centroid
from .utils import computation_wrapper
from scipy.signal import chirp, find_peaks, peak_widths


__all__ = ['gaussian_sigma_width', 'gaussian_fwhm', 'fwhm', 'fwzi']


def gaussian_sigma_width(spectrum, regions=None):
    """
    Estimate the width of the spectrum using a second-moment analysis.

    The value is scaled to match the sigma/standard deviation parameter of a
    standard Gaussian profile. This will be calculated over the regions, if
    they are specified.

    Parameters
    ----------
    spectrum : `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object over which the width will be calculated.

    regions: `~specutils.utils.SpectralRegion` or list of `~specutils.utils.SpectralRegion`
        Region within the spectrum to calculate the gaussian sigma width. If
        regions is `None`, computation is performed over entire spectrum.

    Returns
    -------
    approx_sigma: `~astropy.units.Quantity` or list (based on region input)
        Approximated sigma value of the spectrum

    Notes
    -----
    The spectrum should be continuum subtracted before being passed to this
    function.
    """
    return computation_wrapper(_compute_gaussian_sigma_width, spectrum, regions)


def gaussian_fwhm(spectrum, regions=None):
    """
    Estimate the width of the spectrum using a second-moment analysis.

    The value is scaled to match the full width at half max of a standard
    Gaussian profile.  This will be calculated over the regions, if they are
    specified.

    Parameters
    ----------
    spectrum : `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object overwhich the width will be calculated.

    regions: `~specutils.utils.SpectralRegion` or list of `~specutils.utils.SpectralRegion`
        Region within the spectrum to calculate the FWHM value. If regions is
        `None`, computation is performed over entire spectrum.

    Returns
    -------
    gaussian_fwhm : `~astropy.units.Quantity` or list (based on region input)
        Approximate full width of the signal at half max

    Notes
    -----
    The spectrum should be continuum subtracted before being passed to this
    function.
    """
    return computation_wrapper(_compute_gaussian_fwhm, spectrum, regions)


def fwhm(spectrum, regions=None):
    """
    Compute the true full width half max of the spectrum.

    This makes no assumptions about the shape of the spectrum (e.g. whether it
    is Gaussian). It finds the maximum of the spectrum, and then locates the
    point closest to half max on either side of the maximum, and
    measures the distance between them. This will be calculated over the
    regions, if they are specified.

    Parameters
    ----------
    spectrum : `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object over which the width will be calculated.

    regions: `~specutils.utils.SpectralRegion` or list of `~specutils.utils.SpectralRegion`
        Region within the spectrum to calculate the FWHM value. If regions is
        `None`, computation is performed over entire spectrum.

    Returns
    -------
    whm : `~astropy.units.Quantity` or list (based on region input)
        Full width of the signal at half max

    Notes
    -----
    The spectrum should be continuum subtracted before being passed to this
    function.
    """
    return computation_wrapper(_compute_fwhm, spectrum, regions)


def fwzi(spectrum, regions=None):
    """
    Compute the true full width at zero intensity (i.e. the continuum level)
    of the spectrum.

    This makes no assumptions about the shape of the spectrum. It uses the
    scipy peak-finding to determine the index of the highest flux value, and
    then calculates width at that base of the feature.

    Parameters
    ----------
    spectrum : `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object over which the width will be calculated.

    regions: `~specutils.utils.SpectralRegion` or list of `~specutils.utils.SpectralRegion`
        Region within the spectrum to calculate the FWHM value. If regions is
        `None`, computation is performed over entire spectrum.

    Returns
    -------
    `~astropy.units.Quantity` or list (based on region input)
        Full width of the signal at zero intensity.

    Notes
    -----
    The spectrum must be continuum subtracted before being passed to this
    function.
    """
    return computation_wrapper(_compute_fwzi, spectrum, regions)


def _compute_fwzi(spectrum, regions=None):
    if regions is not None:
        calc_spectrum = extract_region(spectrum, regions)
    else:
        calc_spectrum = spectrum

    # Create a copy of the flux array to ensure the value on the spectrum
    # object is not altered.
    disp, flux = calc_spectrum.spectral_axis, calc_spectrum.flux.copy()

    # For noisy data, ensure that the search from the centroid stops on
    # either side once the flux value reaches zero.
    flux[flux < 0] = 0

    def find_width(data):
        # Find the peaks in the flux data
        peaks, _ = find_peaks(data)

        # Find the index of the maximum peak value in the found peak list
        peak_ind = [peaks[np.argmin(np.abs(
            np.array(peaks) - np.argmin(np.abs(data - np.max(data)))))]]

        # Calculate the width for the given feature
        widths, _, _, _ = \
            peak_widths(data, peak_ind, rel_height=1-1e-7)

        return widths[0] * disp.unit

    if flux.ndim > 1:
        tot_widths = []

        for i in range(flux.shape[0]):
            tot_widths.append(find_width(flux[i]))

        return tot_widths

    return find_width(flux)


def _compute_gaussian_fwhm(spectrum, regions=None):
    """
    This is a helper function for the above `gaussian_fwhm()` method.
    """

    fwhm = _compute_gaussian_sigma_width(spectrum, regions) * gaussian_sigma_to_fwhm

    return fwhm


def _compute_gaussian_sigma_width(spectrum, regions=None):
    """
    This is a helper function for the above `gaussian_sigma_width()` method.
    """

    if regions is not None:
        calc_spectrum = extract_region(spectrum, regions)
    else:
        calc_spectrum = spectrum

    flux = calc_spectrum.flux
    spectral_axis = calc_spectrum.spectral_axis

    centroid_result = centroid(spectrum, regions)

    if flux.ndim > 1:
        spectral_axis = np.broadcast_to(spectral_axis, flux.shape, subok=True)
        centroid_result = centroid_result[:, np.newaxis]

    dx = spectral_axis - centroid_result
    sigma = np.sqrt(np.sum((dx * dx) * flux, axis=-1) / np.sum(flux, axis=-1))

    return sigma


def _compute_single_fwhm(flux, spectral_axis):

    argmax = np.argmax(flux)
    halfval = flux[argmax] / 2

    left = flux[:argmax] <= halfval
    right = flux[argmax+1:] <= halfval

    # Highest signal at the first point
    if np.sum(left) == 0:
        l_idx = 0
    else:
        l_idx = np.where(left == True)[0][-1]

    # Highest signal at the last point
    if np.sum(right) == 0:
        r_idx = len(flux)-1
    else:
        r_idx = np.where(right == True)[0][0] + argmax

    return spectral_axis[r_idx] - spectral_axis[l_idx]


def _compute_fwhm(spectrum, regions=None):
    """
    This is a helper function for the above `fwhm()` method.
    """

    if regions is not None:
        calc_spectrum = extract_region(spectrum, regions)
    else:
        calc_spectrum = spectrum

    flux = calc_spectrum.flux
    spectral_axis = calc_spectrum.spectral_axis

    if flux.ndim > 1:
        return [_compute_single_fwhm(x, spectral_axis) for x in flux]
    else:
        return _compute_single_fwhm(flux, spectral_axis)
