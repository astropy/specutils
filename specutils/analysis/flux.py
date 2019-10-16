"""
A module for analysis tools focused on determining fluxes of spectral features.
"""

import numpy as np
from ..manipulation import extract_region
from .utils import computation_wrapper
from astropy.stats import sigma_clip
from astropy.stats import median_absolute_deviation


__all__ = ['line_flux', 'equivalent_width', 'is_continuum_near_zero']


def line_flux(spectrum, regions=None):
    """
    Computes the integrated flux in a spectrum or region of a spectrum.

    Applies to the whole spectrum by default, but can be limited to a specific
    feature (like a spectral line) if a region is given.

    Parameters
    ----------
    spectrum : Spectrum1D
        The spectrum object over which the summed flux will be calculated.

    regions : `~specutils.SpectralRegion` or list of `~specutils.SpectralRegion`
        Region within the spectrum to calculate the gaussian sigma width. If
        regions is `None`, computation is performed over entire spectrum.

    Returns
    -------
    flux : `~astropy.units.Quantity`
        Flux in the provided spectrum (or regions). Unit isthe ``spectrum``'s'
        ``flux`` unit times ``spectral_axis`` unit.

    Notes
    -----
    While the flux can be computed on any spectrum or region, it should be
    continuum-subtracted to compute actual line fluxes.
    """
    return computation_wrapper(_compute_line_flux, spectrum, regions)


def equivalent_width(spectrum, continuum=1, regions=None):
    """
    Computes the equivalent width of a region of the spectrum.

    Applies to the whole spectrum by default, but can be limited to a specific
    feature (like a spectral line) if a region is given.

    Parameters
    ----------
    spectrum : Spectrum1D
        The spectrum object overwhich the equivalent width will be calculated.

    regions: `~specutils.SpectralRegion` or list of `~specutils.SpectralRegion`
        Region within the spectrum to calculate the gaussian sigma width. If
        regions is `None`, computation is performed over entire spectrum.

    continuum : ``1`` or `~astropy.units.Quantity`, optional
        Value to assume is the continuum level.  For the special value ``1``
        (without units), ``1`` in whatever the units of the ``spectrum.flux``
        will be assumed, otherwise units are required and must be the same as
        the ``spectrum.flux``.

    Returns
    -------
    ew : `~astropy.units.Quantity`
        Equivalent width calculation, in the same units as the ``spectrum``'s
        ``spectral_axis``.

    Notes
    -----
    To do a standard equivalent width measurement, the ``spectrum`` should be
    continuum-normalized to whatever ``continuum`` is before this function is
    called.

    """

    kwargs = dict(continuum=continuum)
    return computation_wrapper(_compute_equivalent_width, spectrum, regions, **kwargs)


def _compute_line_flux(spectrum, regions=None):

    if regions is not None:
        calc_spectrum = extract_region(spectrum, regions)
    else:
        calc_spectrum = spectrum

    # Average dispersion in the line region
    avg_dx = np.diff(calc_spectrum.spectral_axis)

    line_flux = np.sum(calc_spectrum.flux[1:] * avg_dx)

    # TODO: we may want to consider converting to erg / cm^2 / sec by default
    return line_flux


def _compute_equivalent_width(spectrum, continuum=1, regions=None):

    if regions is not None:
        calc_spectrum = extract_region(spectrum, regions)
    else:
        calc_spectrum = spectrum

    if continuum == 1:
        continuum = 1*calc_spectrum.flux.unit

    spectral_axis = calc_spectrum.spectral_axis
    dx = spectral_axis[-1] - spectral_axis[0]


    line_flux = _compute_line_flux(spectrum, regions)

    # Calculate equivalent width
    ew =  dx - (line_flux / continuum)

    return ew.to(calc_spectrum.spectral_axis.unit)

def is_continuum_near_zero(spectrum, eps=0.01):
    """
    Determine if the baseline continuum of the spectrum is near zero.

    The value is scaled to match the sigma/standard deviation parameter of a
    standard Gaussian profile. This will be calculated over the regions, if
    they are specified.

    Parameters
    ----------
    spectrum : `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object over which the width will be calculated.

    eps: float
        The tolerance on the quantification to confirm the continuum is
        near zero.

    Returns
    -------
    is_near_zero: bool
        A True if the continuum of the spectrum is with eps, False otherwise.

    Notes
    -----
    The spectrum should be continuum subtracted before being passed to this
    function.
    """

    # If the eps has units then the assumption is that we want
    # to compare the median of the flux regardless of the
    # existence of the uncertainty.
    if hasattr(eps, 'unit'):
        return np.median(spectrum.flux) < eps

    # If eps does not have a unit, ie it is not a quantity, then
    # we are going to calculate based on the S/N if the uncertainty
    # exists.
    if hasattr(spectrum, 'uncertainty') and spectrum.uncertainty:
        return np.median(spectrum.flux / spectrum.uncertainty.quantity) < eps
    else:
        raise Exception('Spectrum flux has units, eps does not, either include uncertainty or use units on eps')
        #return np.median(spectrum.flux) / median_absolute_deviation(spectrum.flux) < eps
