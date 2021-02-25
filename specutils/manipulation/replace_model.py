import numpy as np

from scipy.interpolate import CubicSpline

from ..spectra import Spectrum1D


def model_replace(spectrum, spline_knots, extrapolation_treatment='data_fill',
                   interpolate_uncertainty=True):
    """
    Generates a new spectrum with a region replaced by a smooth spline.

    Parameters
    ----------

    spectrum : `~specutils.Spectrum1D`
        The spectrum to be modified.

    spline_knots : 1d `~astropy.units.Quantity`
        List of spectral axis values to be used as spline knots.
        They all should share the same units, which can be different
        from the units of the input spectrum spectral axis.

    extrapolation_treatment : str
        What to do with data off the edges of the region encompassed by
        the spline knots. Default is ``'data_fill'`` to have points filled
        with the input flux values. ``'zero_fill'`` sets them to zero.

    interpolate_uncertainty : bool
        If True, the uncertainty is also interpolated over the replaced
        segment.

    Returns
    -------
    spectrum : `~specutils.Spectrum1D`
        The spectrum with the region replaced by spline values, and data
        values or zeros outside the spline region. The spectral axis will
        have the same units as the spline knots.

    Examples
    --------

    >>> import numpy as np
    >>> import astropy.units as u
    >>> from specutils.spectra.spectrum1d import Spectrum1D
    >>> from specutils.manipulation.replace_segment import spline_replace
    >>> wave_val = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    >>> flux_val = np.array([2, 4, 6, 8, 10, 12, 14, 16, 18, 20])
    >>> input_spectrum = Spectrum1D(spectral_axis=wave_val * u.AA, flux=flux_val * u.mJy)
    >>> spline_knots = [3.5, 4.7, 6.8, 7.1] * u.AA
    >>> result = spline_replace(input_spectrum, spline_knots) # doctest: +IGNORE_OUTPUT

    """
    if extrapolation_treatment not in ('data_fill', 'zero_fill'):
        raise ValueError('invalid extrapolation_treatment value: ' + str(extrapolation_treatment))

    # Output spectral axis will have the same units as spline knots
    new_spectral_axis = spectrum.spectral_axis.to(spline_knots.unit)

    # Interpolate on input data to get flux at each spline knot
    flux_spline = CubicSpline(new_spectral_axis.value, spectrum.flux.value,
                              extrapolate=False)
    spline_knots_flux_val = flux_spline(spline_knots.value)

    # Now, evaluate spline over spline knots only, and use it
    # to compute flux at each spectral value in input spectrum.
    flux_spline_2 = CubicSpline(spline_knots.value, spline_knots_flux_val,
                                extrapolate=False)
    out_flux_val = flux_spline_2(new_spectral_axis.value)

    # Do the same in case we want to interpolate the uncertainty.
    # Otherwise, do not propagate uncertainty into output.
    out_uncert_val = None
    if spectrum.uncertainty is not None and interpolate_uncertainty:
        uncert_spline = CubicSpline(new_spectral_axis.value, spectrum.uncertainty.quantity,
                                  extrapolate=False)
        spline_knots_uncert_val = uncert_spline(spline_knots.value)

        uncert_spline_2 = CubicSpline(spline_knots.value, spline_knots_uncert_val,
                                    extrapolate=False)
        out_uncert_val = uncert_spline_2(new_spectral_axis.value)

    # Careful with units handling from here on: astropylts handles the
    # np.where filter differently than the other supported environments.

    # Initialize with zero fill.
    out = np.where(np.isnan(out_flux_val), 0., out_flux_val) * spectrum.flux.unit

    # Fill extrapolated regions with original flux
    if extrapolation_treatment == 'data_fill':
        data = np.where(np.isnan(out_flux_val), spectrum.flux.value, 0.)
        out += (data * spectrum.flux.unit)

    # Do the same in case we want to interpolate the uncertainty.
    # Otherwise, do not propagate uncertainty into output.
    new_unc = None
    if out_uncert_val is not None:
        out_uncert = np.where(np.isnan(out_uncert_val), 0., out_uncert_val) * \
                     spectrum.uncertainty.unit
        data = np.where(np.isnan(out_uncert_val), spectrum.uncertainty.quantity, 0.)
        out_uncert += data

        new_unc = spectrum.uncertainty.__class__(array=out_uncert,
                                                 unit=spectrum.uncertainty.unit)

    return Spectrum1D(spectral_axis=new_spectral_axis, flux=out, uncertainty=new_unc)
