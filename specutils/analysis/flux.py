"""
A module for analysis tools focused on determining fluxes of spectral features.
"""

import numpy as np
from ..manipulation import extract_region
from .utils import computation_wrapper


__all__ = ['line_flux', 'equivalent_width']


def line_flux(spectrum, regions=None):
    """
    Computes the integrated flux in a spectrum or region of a spectrum.

    Applies to the whole spectrum by default, but can be limited to a specific
    feature (such as a spectral line) if a region is given.

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
        Flux in the provided spectrum (or regions). Unit is the ``spectrum``'s'
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
    feature (such as a spectral line) if a region is given.

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

    sign = 1
    # Account for cases with decreasing wavelength
    if calc_spectrum.spectral_axis[0] > calc_spectrum.spectral_axis[-1]:
        sign = -1

    flux = np.trapz(calc_spectrum.flux, x=calc_spectrum.spectral_axis.value)
    flux_units = calc_spectrum.flux.unit * calc_spectrum.spectral_axis.unit
    return flux * sign * flux_units


def _compute_equivalent_width(spectrum, continuum=1, regions=None):

    if regions is not None:
        calc_spectrum = extract_region(spectrum, regions)
    else:
        calc_spectrum = spectrum

    if continuum == 1:
        continuum = 1 * calc_spectrum.flux.unit

    sign = 1
    # Account for cases with decreasing wavelength
    if calc_spectrum.spectral_axis[0] > calc_spectrum.spectral_axis[-1]:
        sign = -1

    ew = np.trapz(1 - calc_spectrum.flux / continuum, x=calc_spectrum.spectral_axis.value)
    ew_units = calc_spectrum.spectral_axis.unit
    return ew * sign * ew_units
