"""
A module for analysis tools focused on determining fluxes of spectral features.
"""

import warnings
from functools import wraps

import numpy as np
from astropy.nddata import VarianceUncertainty

from .. import conf
from ..spectra import Spectrum1D
from ..manipulation import extract_region, LinearInterpolatedResampler
from .utils import computation_wrapper
import astropy.units as u
from astropy.stats import mad_std
from astropy.utils.exceptions import AstropyUserWarning


__all__ = ['line_flux', 'equivalent_width', 'is_continuum_below_threshold',
           'warn_continuum_below_threshold']


def line_flux(spectrum, regions=None,
              mask_interpolation=LinearInterpolatedResampler):
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

    mask_interpolation : ``None`` or `~specutils.manipulation.LinearInterpolatedResampler`
        Interpolator class used to fill up the gaps in the spectrum's flux
        array, when the spectrum mask is not None. If set to ``None``, the
        masked spectral bins are excised from the data without interpolation
        and the bin edges of the adjacent bins are extended to fill the gap.

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
    return computation_wrapper(_compute_line_flux, spectrum, regions,
                               mask_interpolation=mask_interpolation)


def equivalent_width(spectrum, continuum=1, regions=None,
                     mask_interpolation=LinearInterpolatedResampler):
    """
    Computes the equivalent width of a region of the spectrum.

    Applies to the whole spectrum by default, but can be limited to a specific
    feature (like a spectral line) if a region is given.

    Parameters
    ----------
    spectrum : Spectrum1D
        The spectrum object overwhich the equivalent width will be calculated.

    regions: `~specutils.SpectralRegion` or list of `~specutils.SpectralRegion`
        Region within the spectrum to calculate the equivalent width. If
        regions is `None`, computation is performed over entire spectrum.

    continuum : ``1`` or `~astropy.units.Quantity`, optional
        Value to assume is the continuum level.  For the special value ``1``
        (without units), ``1`` in whatever the units of the ``spectrum.flux``
        will be assumed, otherwise units are required and must be the same as
        the ``spectrum.flux``.

    mask_interpolation : ``None`` or `~specutils.manipulation.LinearInterpolatedResampler`
        Interpolator class used to fill up the gaps in the spectrum's flux
        array after an excise operation to ensure the mask shape can always be
        applied when the spectrum mask is not None. If set to ``None``, the
        masked spectral bins are excised from the data without interpolation
        and the bin edges of the adjacent bins are extended to fill the gap.

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
    return computation_wrapper(_compute_equivalent_width, spectrum, regions,
                               mask_interpolation=mask_interpolation, **kwargs)


def _compute_line_flux(spectrum, regions=None,
                       mask_interpolation=LinearInterpolatedResampler):

    if regions is not None:
        calc_spectrum = extract_region(spectrum, regions)
    else:
        calc_spectrum = spectrum

    # Account for the existence of a mask.
    if hasattr(calc_spectrum, 'mask') and calc_spectrum.mask is not None:
        mask = calc_spectrum.mask
        new_spec = Spectrum1D(flux=calc_spectrum.flux[~mask],
                              spectral_axis=calc_spectrum.spectral_axis[~mask])

        if mask_interpolation is None:
            return _compute_line_flux(new_spec)
        else:
            interpolator = mask_interpolation(extrapolation_treatment='zero_fill')
            sp = interpolator(new_spec, calc_spectrum.spectral_axis)
            flux = sp.flux
    else:
        flux = calc_spectrum.flux

    dx = np.abs(np.diff(calc_spectrum.spectral_axis.bin_edges))
    line_flux = np.sum(flux * dx)

    line_flux.uncertainty = None
    line_flux.uncertainty_type = 'stddev'

    if calc_spectrum.uncertainty is not None:
        # Can't handle masks via interpolation here, since interpolators
        # only work with the flux array.
        try:
            variance_q = calc_spectrum.uncertainty.represent_as(VarianceUncertainty).quantity
        except ValueError:
            warnings.warn(f"Uncertainty type "
                          f"'{calc_spectrum.uncertainty.uncertainty_type}' "
                          f"was not recognized by line_flux. Proceeding "
                          f"without uncertainty in result.",
                          AstropyUserWarning)
            variance_q = None

        if variance_q is not None:
            line_flux.uncertainty = np.sqrt(
                np.sum(variance_q * dx**2))

    # TODO: we may want to consider converting to erg / cm^2 / sec by default
    return line_flux


def _compute_equivalent_width(spectrum, continuum=1, regions=None,
                              mask_interpolation=LinearInterpolatedResampler):
    if regions is not None:
        spectrum = extract_region(spectrum, regions)

    # Account for the existence of a mask.
    if hasattr(spectrum, 'mask') and spectrum.mask is not None:
        mask = spectrum.mask
        spectrum = Spectrum1D(
            flux=spectrum.flux[~mask],
            spectral_axis=spectrum.spectral_axis[~mask])

    # Calculate the continuum flux value
    continuum = u.Quantity(continuum, unit=spectrum.flux.unit)

    # If continuum is provided as a single scalar/quantity value, create an
    # array the same size as the processed spectrum. Otherwise, assume that the
    # continuum was provided as an array and use as-is.
    if continuum.size == 1:
        continuum = continuum * np.ones(spectrum.flux.size)

    cont_spec = Spectrum1D(flux=continuum,
                           spectral_axis=spectrum.spectral_axis)
    cont_flux = _compute_line_flux(cont_spec,
                                   mask_interpolation=mask_interpolation)
    line_flux = _compute_line_flux(spectrum,
                                   mask_interpolation=mask_interpolation)

    # Calculate equivalent width
    dx = np.abs(np.diff(spectrum.spectral_axis.bin_edges))
    ew = np.sum((1 - line_flux / cont_flux) * dx)

    # NOTE: uncertainty gets lost during unit translation, so we have to convert first
    ew = ew.to(spectrum.spectral_axis.unit)
    if line_flux.uncertainty is not None:
        # continuum has no uncertainty, here we assume no uncertainty on the bin edges (dx) as well
        # we need to add the line_flux uncertainties in quadrature for the sum above (which was
        # summed over the length of dx)
        ew.uncertainty = np.sqrt((line_flux.uncertainty.value**2)*len(dx)) * dx.unit
    else:
        ew.uncertainty = None

    ew.uncertainty_type = 'stddev'
    return ew


def is_continuum_below_threshold(spectrum, threshold=0.01):
    """
    Determine if the baseline of this spectrum is less than a threshold.
    I.e., an estimate of whether or not the continuum has been subtracted.

    If ``threshold`` is an `~astropy.units.Quantity` with flux units, this
    directly compares the median of the spectrum to the threshold.
    of the flux to the threshold.

    If the threshold is a float or dimensionless quantity then the spectrum's uncertainty will be
    used or an estimate of the uncertainty. If the uncertainty is present then the
    threshold is compared to the median of the flux divided by the
    uncertainty.  If the uncertainty is not present then the threshold
    is compared to the median of the flux divided by the
    `~astropy.stats.mad_std`.

    Parameters
    ----------
    spectrum : `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object over which the width will be calculated.

    threshold: float or `~astropy.units.Quantity`
        The tolerance on the quantification to confirm the continuum is
        near zero.

    Returns
    -------
    is_continuum_below_threshold: bool
        Return True if the continuum of the spectrum is below the threshold, False otherwise.

    """

    flux = spectrum.flux
    uncertainty = spectrum.uncertainty if hasattr(spectrum, 'uncertainty') else None

    # Apply the mask if it exists.
    if hasattr(spectrum, 'mask') and spectrum.mask is not None:
        flux = flux[~spectrum.mask]
        uncertainty = uncertainty[~spectrum.mask] if uncertainty else uncertainty

    # If the threshold has units then the assumption is that we want
    # to compare the median of the flux regardless of the
    # existence of the uncertainty.
    if hasattr(threshold, 'unit') and not threshold.unit == u.dimensionless_unscaled:
        return np.median(flux) < threshold

    # If threshold does not have a unit, ie it is not a quantity, then
    # we are going to calculate based on the S/N if the uncertainty
    # exists.
    if uncertainty and uncertainty.uncertainty_type != 'std':
        return np.median(flux / uncertainty.quantity) < threshold
    else:
        return np.median(flux) / mad_std(flux) < threshold


def warn_continuum_below_threshold(threshold=0.01):
    """
    Decorator for methods that should warn if the baseline
    of the spectrum does not appear to be below a threshold.

    The ``check`` parameter is based on the
    `astropy configuration system <http://docs.astropy.org/en/stable/config/#adding-new-configuration-items>`_.
    Examples are on that page to show how to turn off this type of warning checking.
    """
    def actual_decorator(function):
        @wraps(function)
        def wrapper(*args, **kwargs):
            if conf.do_continuum_function_check:
                spectrum = args[0]
                if not is_continuum_below_threshold(spectrum, threshold):
                    if hasattr(threshold, 'unit'):
                        levelorsnr = 'value'
                    else:
                        levelorsnr = 'signal-to-noise'

                    message = "Spectrum is not below the threshold {} {}. ".format(levelorsnr, threshold)
                    message += "This may indicate you have not continuum subtracted this spectrum (or that you have but it has high SNR features).\n\n"
                    message += ("""If you want to suppress this warning either type """
                                """'specutils.conf.do_continuum_function_check = False' or """
                                """see http://docs.astropy.org/en/stable/config/#adding-new-configuration-items """
                                """for other ways to configure the warning.""")

                    warnings.warn(message, AstropyUserWarning)
            result = function(*args, **kwargs)
            return result
        return wrapper
    return actual_decorator
