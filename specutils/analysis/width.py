"""
A module for analysis tools focused on determining the width of
spectral features.
"""

from astropy.nddata import StdDevUncertainty
from astropy.stats.funcs import gaussian_sigma_to_fwhm
import astropy.uncertainty as unc
import numpy as np
from scipy.signal import find_peaks, peak_widths

from ..manipulation import extract_region
from . import centroid
from .utils import computation_wrapper


__all__ = ['gaussian_sigma_width', 'gaussian_fwhm', 'fwhm', 'fwzi']


def gaussian_sigma_width(spectrum, regions=None, analytic=True):
    """
    Estimate the width of the spectrum using a second-moment analysis.

    The value is scaled to match the sigma/standard deviation parameter of a
    standard Gaussian profile. This will be calculated over the regions, if
    they are specified.

    Parameters
    ----------
    spectrum : `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object over which the width will be calculated.

    regions: `~specutils.SpectralRegion` or list of `~specutils.SpectralRegion`
        Region within the spectrum to calculate the gaussian sigma width. If
        regions is `None`, computation is performed over entire spectrum.

    analytic : bool, optional
        Set this flag to ``False`` to use the `~astropy.uncertainty` distribution-based
        calculation for the width and its uncertainty instead of the default analytic solution.

    Returns
    -------
    approx_sigma: `~astropy.units.Quantity` or list (based on region input)
        Approximated sigma value of the spectrum

    Notes
    -----
    The spectrum should be continuum subtracted before being passed to this
    function.
    """
    return computation_wrapper(_compute_gaussian_sigma_width, spectrum, regions,
                               analytic=analytic)


def gaussian_fwhm(spectrum, regions=None, analytic=True):
    """
    Estimate the width of the spectrum using a second-moment analysis.

    The value is scaled to match the full width at half max of a standard
    Gaussian profile.  This will be calculated over the regions, if they are
    specified.

    Parameters
    ----------
    spectrum : `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object overwhich the width will be calculated.

    regions : `~specutils.SpectralRegion` or list of `~specutils.SpectralRegion`
        Region within the spectrum to calculate the FWHM value. If regions is
        `None`, computation is performed over entire spectrum.

    analytic : bool, optional
        Set this flag to ``False`` to use the `~astropy.uncertainty` distribution-based
        calculation for the fwhm and its uncertainty instead of the default analytic solution.

    Returns
    -------
    gaussian_fwhm : `~astropy.units.Quantity` or list (based on region input)
        Approximate full width of the signal at half max

    Notes
    -----
    The spectrum should be continuum subtracted before being passed to this
    function.
    """
    return computation_wrapper(_compute_gaussian_fwhm, spectrum, regions,
                               analytic=analytic)


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

    regions: `~specutils.SpectralRegion` or list of `~specutils.SpectralRegion`
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

    regions: `~specutils.SpectralRegion` or list of `~specutils.SpectralRegion`
        Region within the spectrum to calculate the FWZI value. If regions is
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
    if hasattr(spectrum, 'mask') and spectrum.mask is not None:
        disp = calc_spectrum.spectral_axis[~spectrum.mask]
        flux = calc_spectrum.flux[~spectrum.mask].copy()
    else:
        disp = calc_spectrum.spectral_axis
        flux = calc_spectrum.flux.copy()

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


def _compute_gaussian_fwhm(spectrum, regions=None, analytic=False):
    """
    This is a helper function for the above `gaussian_fwhm()` method.
    """

    sigma = _compute_gaussian_sigma_width(spectrum, regions, analytic=analytic)

    fwhm = sigma * gaussian_sigma_to_fwhm
    if sigma.uncertainty is not None:
        fwhm.uncertainty = sigma.uncertainty * gaussian_sigma_to_fwhm
    else:
        fwhm.uncertainty = None

    fwhm.uncertainty_type = 'std'
    return fwhm


def _compute_gaussian_sigma_width(spectrum, regions=None, analytic=True):
    """
    This is a helper function for the above `gaussian_sigma_width()` method.
    """

    if not analytic and spectrum.uncertainty is None:
        raise ValueError("Distribution-based calculation can only be used if"
                         " spectrum.uncertainty is not None")

    if regions is not None:
        calc_spectrum = extract_region(spectrum, regions)
    else:
        calc_spectrum = spectrum

    if spectrum.uncertainty is not None:
        flux_uncert = calc_spectrum.uncertainty.represent_as(StdDevUncertainty).quantity
    else:
        # dummy value for uncertainties to avoid extra if-statements when applying mask
        flux_uncert = np.zeros_like(calc_spectrum.flux)

    if hasattr(spectrum, 'mask') and spectrum.mask is not None:
        flux = calc_spectrum.flux[~spectrum.mask]
        spectral_axis = calc_spectrum.spectral_axis[~spectrum.mask].quantity
        flux_uncert = flux_uncert[~calc_spectrum.mask]
    else:
        flux = calc_spectrum.flux
        spectral_axis = calc_spectrum.spectral_axis.quantity

    centroid_result = centroid(spectrum, regions=regions, analytic=analytic)

    if flux.ndim > 1:
        spectral_axis = np.broadcast_to(spectral_axis, flux.shape, subok=True)
        if centroid_result.uncertainty is None:
            temp_uncertainty = None
        else:
            temp_uncertainty = centroid_result.uncertainty[:, np.newaxis]
        centroid_result = centroid_result[:, np.newaxis]
        centroid_result.uncertainty = temp_uncertainty

    if not analytic:
        # Convert the centroid and flux values to astropy uncertainty distributions
        centroid_result = unc.normal(centroid_result, std=centroid_result.uncertainty,
                                     n_samples=1000)
        flux = unc.normal(flux, std=flux_uncert, n_samples=1000)

    dx = (spectral_axis - centroid_result)
    if analytic:
        numerator = np.sum((dx * dx) * flux, axis=-1)
        denom = np.sum(flux, axis=-1)
    else:
        # Account for trailing n_samples axis
        numerator = np.sum((dx * dx) * flux, axis=-2)
        denom = np.sum(flux, axis=-2)

    sigma2 = numerator / denom

    if analytic:
        sigma = np.sqrt(sigma2)
        if centroid_result.uncertainty is not None:
            # NOTE: until/unless disp_uncert is supported, dx_uncert == centroid_result.uncertainty
            disp_uncert = 0.0 * spectral_axis.unit
            dx_uncert = np.sqrt(disp_uncert**2 + centroid_result.uncertainty**2)

            # Uncertainty for each term in the numerator sum
            num_term_uncerts = dx * dx * flux * np.sqrt(2*(dx_uncert/dx)**2 + (flux_uncert/flux)**2)
            # uncertainty (squared) for the numerator, added in quadrature
            num_uncertsq = np.sum(num_term_uncerts**2, axis=-1)
            # uncertainty (squared) for the denomenator
            denom_uncertsq = np.sum(flux_uncert**2)

            sigma2_uncert = numerator/denom * np.sqrt(num_uncertsq * numerator**-2 +
                                                      denom_uncertsq * denom**-2)

            sigma.uncertainty = 0.5 * sigma2_uncert / sigma2 * sigma
        else:
            sigma.uncertainty = None
    else:
        # Need to do pdf_mean and pdf_std here to avoid negatives in sqrt
        sigma2_uncertainty = sigma2.pdf_std()
        sigma2 = sigma2.pdf_mean()
        sigma = np.sqrt(sigma2)
        sigma.uncertainty = 0.5 * sigma * sigma2_uncertainty / sigma2

    sigma.uncertainty_type = 'stddev'
    return sigma


def _compute_single_fwhm(flux, spectral_axis):
    """
    This is a helper function for the above `fwhm()` method.
    """

    # The .value attribute is used here as the following algorithm does not
    # use any array operations and would otherwise introduce a relatively
    # significant overhead factor.  Two-point linear interpolation is used to
    # achieve sub-pixel precision.
    flux_value = flux.value
    spectral_value = spectral_axis.value

    argmax = flux_value.argmax()
    halfval = flux_value[argmax] / 2
    left = flux_value[:argmax] < halfval
    right = flux_value[argmax + 1:] < halfval

    # Highest signal at the first point
    i0 = np.nonzero(left)[0]
    if i0.size == 0:
        left_value = spectral_value[0]
    else:
        i0 = i0[-1]
        i1 = i0 + 1
        left_flux = flux_value[i0]
        left_spectral = spectral_value[i0]
        left_value = ((halfval - left_flux)
                      * (spectral_value[i1] - left_spectral)
                      / (flux_value[i1] - left_flux)
                      + left_spectral)

    # Highest signal at the last point
    i1 = np.nonzero(right)[0]
    if i1.size == 0:
        right_value = spectral_value[-1]
    else:
        i1 = i1[0] + argmax + 1
        i0 = i1 - 1
        left_flux = flux_value[i0]
        left_spectral = spectral_value[i0]
        right_value = ((halfval - left_flux)
                       * (spectral_value[i1] - left_spectral)
                       / (flux_value[i1] - left_flux)
                       + left_spectral)

    return spectral_axis.unit * np.abs(right_value - left_value)


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
