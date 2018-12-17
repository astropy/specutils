from __future__ import division

import numpy as np
import astropy.units as u
from astropy import convolution
from astropy.modeling import models, fitting
from scipy.signal import medfilt

from ..spectra import Spectrum1D

from scipy.interpolate import interp1d
from scipy.optimize import minimize


__all__ = ['align_spectra']


def _interpolate_spectrum(wavelengths, spectral_flux, out_wavelengths):
    """
    Re-interpolate the spectral values onto a new set of
    wavelengths.

    Parameters
    ----------
    wavelengths : `~astropy.units.Quantity` array or numpy array
        The spectral wavelength values of the spectrum

    spectral_flux : `~astropy.units.Quantity` array or numpy array
        The spectral flux values of the spectrum

    out_wavelengths : `~astropy.units.Quantity` array or numpy array
        The spectral wavelengths of the desired output spectrum

    Returns
    -------
    out_spectral_flux
        The interpolated/extrapolated spectral flux values at the
        ``out_wavelengths`` positions.

    """

    interp_func = interp1d(wavelengths, spectral_flux, fill_value='extrapolate')

    return interp_func(out_wavelengths)


def _fit_metric(poly_coefs, old_wave, spec2, out_wave, spec1_scaled):
    """
    Method whose output should be minimized. Used by the fitting
    program in the alignment method.
    """

    #
    # Evaluate the coefficients
    #

    altered_wave2 = np.polyval(poly_coefs, old_wave)

    #
    # Interpolate the spectrum at the new altered wavelengths.
    #

    altered_spec2 = _interpolate_spectrum(altered_wave2, spec2, out_wave) #resample the spectrum

    #
    # Calculate the difference between the altered spectrum
    # and the scaled.
    #

    diff_spec = altered_spec2 - spec1_scaled

    return diff_spec.std()


def align_spectra(spectrum1, spectrum2):
    """
    Align the features in ``spectrum2`` to the features in ``spectrum1``.

    Parameters
    ----------
    spectrum1 : `~specutils.Spectrum1D`
        Base spectrum.

    spectrum2 : `~specutils.Spectrum1D`
        Spectrum to align to the base spectrum.

    Returns
    -------
    out_spectrum: `~specutils.Spectrum1D`
        The interpolated/extrapolated spectral flux values at the
        ``out_wavelengths`` positions.

    """

    #
    # Grab the underlying data.
    #

    wave1, spec1 = spectrum1.spectral_axis, spectrum1.flux
    wave2, spec2 = spectrum2.spectral_axis, spectrum2.flux

    #
    # Since the spectra might be of different shapes, we need to isolate the portions
    # of the spectra where they overlap
    #

    wave_lo, wave_hi = max(wave1.min(), wave2.min()), min(wave1.max(), wave2.max())
    wave1_overlap = (wave1 >= wave_lo) & (wave1 <= wave_hi) & (spec1 != 0.0)
    wave2_overlap = (wave2 >= wave_lo) & (wave2 <= wave_hi) & (spec2 != 0.0)

    #
    # We also need to normalize the fluxes to make our wavelength fits
    # better.
    #

    out_wave = wave1[wave1_overlap]
    interp_spec2 = _interpolate_spectrum(wave2[wave2_overlap], spec2[wave2_overlap], out_wave)
    normalized_spec = interp_spec2 / spec1[wave1_overlap]

    #
    # If the spectra are not already well-aligned, more robust fitting
    # methods, probably including outlier rejection, would be necessary.
    #

    norm_model = models.Polynomial1D(2)
    fitter = fitting.LinearLSQFitter()

    ok_px = (np.isfinite(normalized_spec))
    normalization = fitter(norm_model, out_wave[ok_px], normalized_spec[ok_px])

    flux_scale_factor = normalization(out_wave)
    spec1_scaled = spec1[wave1_overlap] * flux_scale_factor

    #
    # Try a 2nd-degree polynomial to start; the coefficient array can
    # be any size, depending on what degree polynomial you wish to fit.
    #

    pixel_fit = minimize(_fit_metric, np.array([0., 1., 0.]),
                         args=(wave2.value[wave2_overlap],spec2[wave2_overlap], out_wave, spec1_scaled, ), method='Nelder-Mead')

    alt_wave2 = np.polyval(pixel_fit.x, wave2.value[wave2_overlap])
    alt_spec2 = _interpolate_spectrum(alt_wave2, spec2[wave2_overlap], out_wave)

    #
    # Return a new object with the smoothed flux.
    #

    return Spectrum1D(flux=u.Quantity(alt_spec2, spectrum1.unit),
                      spectral_axis=u.Quantity(alt_wave2, spectrum1.spectral_axis_unit))
