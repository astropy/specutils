from scipy.interpolate import CubicSpline

from ..spectra import Spectrum1D


def spline_replace(spectrum, spline_points, extrapolation_treatment='nan_fill'):

    if extrapolation_treatment not in ('nan_fill', 'zero_fill'):
        raise ValueError('invalid extrapolation_treatment value: ' + str(extrapolation_treatment))

    # Output spectral axis will have the same units as spline points
    new_spectral_axis = spectrum.spectral_axis.to(spline_points.unit)

    # Interpolate on input data to get flux at each spline point
    flux_spline = CubicSpline(new_spectral_axis.value, spectrum.flux.value,
                              extrapolate=extrapolation_treatment != 'nan_fill')
    spline_points_flux_val = flux_spline(spline_points.value)

    # Now, evaluate spline over spline points only, and use it
    # to compute flux at each spectral value in input spectrum.
    flux_spline_2 = CubicSpline(spline_points.value, spline_points_flux_val,
                                extrapolate=extrapolation_treatment != 'nan_fill')
    out_flux_val = flux_spline_2(new_spectral_axis.value) * spectrum.flux.unit

    return Spectrum1D(spectral_axis=new_spectral_axis, flux=out_flux_val)
