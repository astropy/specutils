import numpy as np
from astropy.stats.funcs import gaussian_fwhm_to_sigma
from ..spectra import SpectralRegion


__all__ = ['gaussian_sigma_width', 'gaussian_fwhm', 'fwhm']


def _computation_wrapper(func, spectrum, region):

    # No region, therefore whole spectrum.
    if region is None:
        return func(spectrum)

    # Single region
    elif isinstance(region, SpectralRegion):
        return func(spectrum, region=region)

    # List of regions
    elif isinstance(region, list):
        return [func(spectrum, region=reg) for reg in region]


def gaussian_sigma_width(spectrum, region=None):
    """
    Estimate the full width of the spectrum based on an approximate value of
    sigma.  This will be calculated over the regions, if they are specified.

    Parameters
    ----------
    spectrum : `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object over which the width will be calculated.

    region: `~specutils.utils.SpectralRegion` or list of `~specutils.utils.SpectralRegion`
        Region within the spectrum to calculate the gaussian sigma width.

    Returns
    -------
    full_width : `~astropy.units.Quantity` or list (based on region input)
        Approximate full width of the signal
    """
    return _computation_wrapper(_compute_gaussian_sigma_width, spectrum, region)


def gaussian_fwhm(spectrum, region=None):
    """
    Estimate the full width half max of the spectrum assuming it is gaussian.
    This will be calculated over the regions, if they are specified.

    Parameters
    ----------
    spectrum : `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object overwhich the width will be calculated.

    region: `~specutils.utils.SpectralRegion` or list of `~specutils.utils.SpectralRegion`
        Region within the spectrum to calculate the FWHM value.

    Returns
    -------
    gaussian_fwhm : `~astropy.units.Quantity` or list (based on region input)
        Approximate full width of the signal at half max
    """
    return _computation_wrapper(_compute_gaussian_fwhm, spectrum, region)


def fwhm(spectrum, region=None):
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

    region: `~specutils.utils.SpectralRegion` or list of `~specutils.utils.SpectralRegion`
        Region within the spectrum to calculate the FWHM value.

    Returns
    -------
    whm : `~astropy.units.Quantity` or list (based on region input)
        Full width of the signal at half max
    """
    return _computation_wrapper(_compute_fwhm, spectrum, region)


def _compute_gaussian_fwhm(spectrum, region=None):
    """
    Estimate the full width of the spectrum at half max.

    Parameters
    ----------
    spectrum : `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object overwhich the width will be calculated.

    region: `~specutils.utils.SpectralRegion`
        Region within the spectrum to calculate the FWHM value.

    Returns
    -------
    gaussian_fwhm : `~astropy.units.Quantity` or list (based on region input)
        Approximate full width of the signal at half max

    Notes
    -----
    This is a helper function for the above `gaussian_fwhm()` method.

    """

    if region is not None:
        calc_spectrum = region.extract(spectrum)
    else:
        calc_spectrum = spectrum

    flux = calc_spectrum.flux
    spectrum = calc_spectrum.spectral_axis

    dx = spectrum - np.mean(spectrum)
    fwhm = 2 * np.sqrt(np.sum((dx * dx) * flux, axis=-1) / np.sum(flux, axis=-1))

    return fwhm


def _compute_gaussian_sigma_width(spectrum, region=None):
    """
    Estimate the full width of the spectrum based on an approximate value of
    sigma.

    Parameters
    ----------
    spectrum : `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object overwhich the width will be calculated.

    region: `~specutils.utils.SpectralRegion`
        Region within the spectrum to calculate the gaussian sigma width.

    Returns
    -------
    full_width : `~astropy.units.Quantity` or list (based on region input)
        Approximate full width of the signal

    Notes
    -----
    This is a helper function for the above `gaussian_sigma_width()` method.

    """

    fwhm = _compute_gaussian_fwhm(spectrum, region)
    sigma = fwhm * gaussian_fwhm_to_sigma

    return sigma * 2


def _arghalf(array, halfval):
    """
    Find index of the item closest to the half max
    """
    return np.abs(array - halfval).argmin()


def _compute_single_fwhm(flux, spectrum):

    argmax = np.argmax(flux)
    halfval = flux[argmax] / 2

    max_idx = len(flux) - 1
    l_idx = _arghalf(flux[:argmax], halfval) if argmax > 0 else 0
    r_idx = _arghalf(flux[argmax+1:], halfval) + argmax+1 if argmax < max_idx else max_idx

    return spectrum[r_idx] - spectrum[l_idx]


def _compute_fwhm(spectrum, region=None):
    """
    Calculate the full width of the spectrum at half max.

    Parameters
    ----------
    spectrum : `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object over which the width will be calculated.

    region: `~specutils.utils.SpectralRegion`
        Region within the spectrum to calculate FWHMh.

    Returns
    -------
    fwhm : `~astropy.units.Quantity` or list (based on region input)
        Computed full width of the signal at half max

    Notes
    -----
    This is a helper function for the above `fwhm()` method.

    """

    if region is not None:
        calc_spectrum = region.extract(spectrum)
    else:
        calc_spectrum = spectrum

    flux = calc_spectrum.flux
    spectrum = calc_spectrum.spectral_axis

    if flux.ndim > 1:
        return [_compute_single_fwhm(x, spectrum) for x in flux]
    else:
        return _compute_single_fwhm(flux, spectrum)
