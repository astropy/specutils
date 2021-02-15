import numpy as np

from scipy.interpolate import CubicSpline

from ..spectra import Spectrum1D


def spline_replace(spectrum, spline_knots, extrapolation_treatment='zero_fill'):
    """
    Generates a new spectrum with a region replaced by a smooth spline.

    Parameters
    ----------

    spectrum : `~specutils.Spectrum1D`
        The spectrum to be modified.

    spline_knots : array of `~astropy.units.Quantity`
        List of spectral axis values to be used as spline knots.
        They all should share the same units, which can be different
        from the units of the input spectrum spectral axis.

    extrapolation_treatment : str
        What to do with data off the edges of the region encompassed by
        the spline knots. Can be ``'data_fill'`` to have points filled
        with the input flux values, or ``'zero_fill'`` to be set to zero.

    Returns
    -------
    spectrum : `~specutils.Spectrum1D`
        The spectrum with the region replaced by spline values, and data
        values or zeros outside the spline region. The spectral axis will
        have the same units as the spline knots.

    Notes
    -----
    The uncertainty attribute, if present in the input spectrum, is not
    carried over into the output.

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
    >>> result = spline_replace(input_spectrum, spline_knots,
    ...     extrapolation_treatment='data_fill') # doctest: +IGNORE_OUTPUT

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
    out_flux_val = flux_spline_2(new_spectral_axis.value) * spectrum.flux.unit

    # default behavior (zero_fill)
    out = np.where(np.isnan(out_flux_val), 0., out_flux_val)

    if extrapolation_treatment == 'data_fill':
        data = np.where(np.isnan(out_flux_val), spectrum.flux.value, 0.)
        out += data * spectrum.flux.unit

    return Spectrum1D(spectral_axis=new_spectral_axis, flux=out)
