import numpy as np

from scipy.interpolate import CubicSpline
from astropy.units import Quantity
from astropy.modeling import Fittable1DModel

from ..spectra import Spectrum1D
from ..utils import QuantityModel
from . import extract_region


def model_replace(spectrum, replace_region, model=10, extrapolation_treatment='data_fill',
                  interpolate_uncertainty=True):
    """
    Generates a new spectrum with a region replaced by a smooth spline.

    Parameters
    ----------

    spectrum : `~specutils.Spectrum1D`
        The spectrum to be modified.

    replace_region : `~specutils.SpectralRegion`
        Spectral region that specifies the region to be replaced. If None,
        parameter `model` below should define spline knots explicitly in the
        form of an `~astropy.units.Quantity` list or array

    model :
        An `~astropy.modeling` model object, which is assumed to have
        already been fit to the spectrum;
        Or
        a list or array of spectral axis values to be used as spline knots.
        They all should share the same units, which can be different
        from the units of the input spectrum spectral axis, but most
        be of compatible physical type;
        Or
        an integer value that will be used to build a list of equally-spaced
        knots, based on the `replace_region` instance.

    extrapolation_treatment : str
        What to do with data off the edges of the region encompassed by
        the spline knots. Default is ``'data_fill'`` to have points filled
        with the input flux values. ``'zero_fill'`` sets them to zero.

    interpolate_uncertainty : bool
        If True, the uncertainty, if present in the input spectrum, is
        also interpolated over the replaced segment.

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
    >>> from specutils.manipulation.model_replace import model_replace
    >>> wave_val = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    >>> flux_val = np.array([2, 4, 6, 8, 10, 12, 14, 16, 18, 20])
    >>> input_spectrum = Spectrum1D(spectral_axis=wave_val * u.AA, flux=flux_val * u.mJy)
    >>> spline_knots = [3.5, 4.7, 6.8, 7.1] * u.AA
    >>> model_replace(input_spectrum, None, spline_knots)
    <Spectrum1D(flux=<Quantity [ 2.,  4.,  6.,  8., 10., 12., 14., 16., 18., 20.] mJy> (shape=(10,), mean=11.00000 mJy); spectral_axis=<SpectralAxis [ 1.  2.  3. ...  8.  9. 10.] Angstrom> (length=10))>

    """
    if extrapolation_treatment not in ('data_fill', 'zero_fill'):
        raise ValueError('invalid extrapolation_treatment value: ' +
                         str(extrapolation_treatment))

    # If input model is an array with spline knots and region
    # is not present, use knots directly to fit spline.
    if isinstance(model, Quantity) and replace_region is None:
        new_spectral_axis = spectrum.spectral_axis.to(model.unit)

        out_flux_val, out_uncert_val = _compute_spline_values(spectrum,
                                                              model,
                                                              new_spectral_axis,
                                                              interpolate_uncertainty)

    # If input model is an int, use it and the spectral region to build a
    # list of equally spaced spline knots.
    elif isinstance(model, int) and replace_region is not None:
        nknots = max(model, 3)
        dw = (replace_region.upper - replace_region.lower) / (nknots - 1)

        spline_knots = []
        for k in range(nknots):
            spline_knots.append(replace_region.lower + k * dw)

        spline_knots = spline_knots * replace_region.lower.unit
        new_spectral_axis = spectrum.spectral_axis.to(spline_knots.unit)

        out_flux_val, out_uncert_val = _compute_spline_values(spectrum,
                                                              spline_knots,
                                                              new_spectral_axis,
                                                              interpolate_uncertainty)

    # If input model is a fitted model, use it and the spectral region to place the
    # model values over the relevant stretch of the spectrum's spectral axis.
    elif (isinstance(model, Fittable1DModel) or (isinstance(model, QuantityModel))) \
            and replace_region is not None:
        new_spectral_axis = spectrum.spectral_axis

        subspectrum = extract_region(spectrum, replace_region)

        model_values = model(subspectrum.spectral_axis)

        out_flux_val = np.full(spectrum.spectral_axis.shape, np.nan)
        indices_up = np.where(spectrum.spectral_axis >= replace_region.lower)
        indices_dn = np.where(spectrum.spectral_axis <= replace_region.upper)
        i = int(indices_up[0][0])
        j = int(indices_dn[0][-1]) + 1

        out_flux_val[i:j] = model_values

        # models do not propagate uncertainties.
        out_uncert_val = None

    else:
        raise NotImplementedError("This combination of input parameters is not yet implemented.")

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
        data = np.where(np.isnan(out_uncert_val), spectrum.uncertainty.quantity.value, 0.)
        out_uncert += (data * spectrum.uncertainty.unit)

        new_unc = spectrum.uncertainty.__class__(array=out_uncert,
                                                 unit=spectrum.uncertainty.unit)

    return Spectrum1D(spectral_axis=new_spectral_axis, flux=out, uncertainty=new_unc)


def _compute_spline_values(spectrum, spline_knots, new_spectral_axis, interpolate_uncertainty):

    # Compute output flux values interpolated over the spline knots.
    out_flux_val = _interpolate_spline(spectrum.flux.value, new_spectral_axis, spline_knots)

    # Do the same in case we want to interpolate the uncertainty.
    # Otherwise, do not propagate uncertainty into output.
    out_uncert_val = None
    if spectrum.uncertainty is not None and interpolate_uncertainty:
        out_uncert_val = _interpolate_spline(spectrum.uncertainty.quantity,
                                             new_spectral_axis, spline_knots)

    return out_flux_val, out_uncert_val


def _interpolate_spline(input_values, spectral_axis, spline_knots):

    # Create spline to interpolate on input data.
    # Knots are the spectral axis values themselves.
    spline_1 = CubicSpline(spectral_axis.value, input_values,
                           extrapolate=False)

    # Now use that spline interpolator to compute interpolated
    # values from the input array at each spline knot.
    values_at_spline_knots = spline_1(spline_knots.value)

    # Finally, compute another spline interpolator over only the values
    # at spline knots, and use it to interpolate the output value at each
    # point on the spectral axis.
    spline_2 = CubicSpline(spline_knots.value, values_at_spline_knots,
                           extrapolate=False)

    return spline_2(spectral_axis.value)
