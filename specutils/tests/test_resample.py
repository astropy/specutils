import numpy as np
import pytest
import astropy.units as u
from astropy.nddata import VarianceUncertainty, InverseVariance, StdDevUncertainty
from astropy.tests.helper import assert_quantity_allclose

from ..spectra.spectrum1d import Spectrum1D, SpectralAxis
from ..manipulation.resample import FluxConservingResampler, LinearInterpolatedResampler, SplineInterpolatedResampler


@pytest.fixture(params=[FluxConservingResampler, LinearInterpolatedResampler, SplineInterpolatedResampler])
def all_resamplers(request):
    return request.param

# todo: Should add tests for different weighting options once those
# are more solidified.


def test_same_grid_fluxconserving(simulated_spectra):
    """
    Test that feeding in the original dispersion axis returns the
    same flux after resampling.
    """
    input_spectra = simulated_spectra.s1_um_mJy_e1
    input_spectra.uncertainty = InverseVariance([0.5]*len(simulated_spectra.s1_um_mJy_e1.flux))

    inst = FluxConservingResampler()
    results = inst(input_spectra,
                   simulated_spectra.s1_um_mJy_e1.spectral_axis)

    assert np.allclose(np.array(simulated_spectra.s1_um_mJy_e1.flux),
                       np.array(results.flux))
    assert np.allclose(input_spectra.uncertainty.array,
                       results.uncertainty.array)


def test_expanded_grid_fluxconserving():
    """
    New dispersion axis has more bins then input dispersion axis
    """

    flux_val = np.array([1, 3, 7, 6, 20])
    wave_val = np.array([2, 4, 12, 16, 20])
    input_spectra = Spectrum1D(flux=flux_val * u.mJy, spectral_axis=wave_val * u.nm)
    resamp_grid = [1, 5, 9, 13, 14, 17, 21, 22, 23] * u.nm

    inst = FluxConservingResampler()
    results = inst(input_spectra, resamp_grid)

    assert_quantity_allclose(results.flux,
                             np.array([1., 3., 6., 7., 6.25, 10., 20.,
                                       np.nan, np.nan])*u.mJy)


def test_stddev_uncert_propogation():
    """
    Check uncertainty propagation if input uncertainty is InverseVariance
    """
    flux_val = np.array([1, 3, 7, 6, 20])
    wave_val = np.array([20, 30, 40, 50, 60])
    input_spectra = Spectrum1D(flux=flux_val * u.mJy, spectral_axis=wave_val * u.AA,
                               uncertainty=StdDevUncertainty([0.1, 0.25, 0.1, 0.25, 0.1]))

    inst = FluxConservingResampler()
    results = inst(input_spectra, [25, 35, 50, 55]*u.AA)

    assert np.allclose(results.uncertainty.array,
                       np.array([31.59810022, 38.89549208, 21.30305717, 31.59810022]))


def delta_wl(saxis):
    """
    A helper function that computes the "size" of a bin given the bin centers
    for testing the flux conservation
    """
    l_widths = (saxis[1] - saxis[0])
    r_widths = (saxis[-1] - saxis[-2])
    # if three bins 0,1,2; want width of central bin.  width is avg of 1/2 minus
    # average of 0/1: (i1 + i2)/2 - (i0 + i1)/2 = (i2 - i0)/2
    mid_widths = (saxis[2:] - saxis[:-2]) / 2

    return np.concatenate([[l_widths.value], mid_widths.value, [r_widths.value]])*saxis.unit


@pytest.mark.parametrize("specflux,specwavebins,outwavebins", [
    ([1, 3, 2], [4000, 5000, 6000, 7000], np.linspace(4000, 7000, 5)),
    ([1, 3, 2, 1], np.linspace(4000, 7000, 5), [4000, 5000, 6000, 7000])])
def test_flux_conservation(specflux, specwavebins, outwavebins):
    """
    A few simple cases to programatically ensure flux is conserved in the
    resampling algorithm
    """
    specwavebins = specwavebins*u.AA
    outwavebins = outwavebins*u.AA
    specflux = specflux*u.AB

    specwave = (specwavebins[:-1] + specwavebins[1:])/2
    outwave = (outwavebins[:-1] + outwavebins[1:])/2

    in_spec = Spectrum1D(spectral_axis=specwave, flux=specflux)
    out_spec = FluxConservingResampler()(in_spec, outwave)

    in_dwl = delta_wl(in_spec.spectral_axis)
    out_dwl = delta_wl(out_spec.spectral_axis)

    flux_in = np.sum(in_spec.flux * in_dwl)
    flux_out = np.sum(out_spec.flux * out_dwl)

    assert_quantity_allclose(flux_in, flux_out)


def test_multi_dim_spectrum1D():
    """
    Test for input spectrum1Ds that have a two dimensional flux and
    uncertainty.
    """
    flux_2d = np.array([np.ones(10) * 5, np.ones(10) * 6, np.ones(10) * 7])

    input_spectra = Spectrum1D(spectral_axis=np.arange(5000, 5010) * u.AA,
                               flux=flux_2d * u.Jy,
                               uncertainty=StdDevUncertainty(flux_2d / 10))

    inst = FluxConservingResampler()
    results = inst(input_spectra, [5001, 5003, 5005, 5007] * u.AA)

    assert_quantity_allclose(results.flux,
                             np.array([[5., 5., 5., 5.],
                                       [6., 6., 6., 6.],
                                       [7., 7., 7., 7.]]) * u.Jy)
    assert np.allclose(results.uncertainty.array,
                       np.array([[6.53197265, 6.53197265, 6.53197265, 6.53197265],
                                 [4.53609212, 4.53609212, 4.53609212, 4.53609212],
                                 [3.33263911, 3.33263911, 3.33263911, 3.33263911]]))


def test_expanded_grid_interp_linear():
    """
    New dispersion axis has more bins then input dispersion axis
    """

    flux_val = np.array([1, 3, 7, 6, 20])
    wave_val = np.array([2, 4, 12, 16, 20])
    input_spectra = Spectrum1D(spectral_axis=wave_val * u.AA, flux=flux_val * u.mJy)
    resamp_grid = [1, 5, 9, 13, 14, 17, 21, 22, 23] * u.AA

    inst = LinearInterpolatedResampler()
    results = inst(input_spectra, resamp_grid)

    assert_quantity_allclose(results.flux,
                             np.array([np.nan, 3.5, 5.5, 6.75, 6.5, 9.5, np.nan,
                                       np.nan, np.nan])*u.mJy)

    inst = LinearInterpolatedResampler(extrapolation_treatment="truncate")
    results = inst(input_spectra, resamp_grid)

    assert_quantity_allclose(results.flux,
                             np.array([3.5, 5.5, 6.75, 6.5, 9.5])*u.mJy)


def test_expanded_grid_interp_spline():
    """
    New dispersion axis has more bins then input dispersion axis
    """

    flux_val = np.array([1, 3, 7, 6, 20])
    wave_val = np.array([2, 4, 12, 16, 20])
    input_spectra = Spectrum1D(spectral_axis=wave_val * u.AA, flux=flux_val * u.mJy)
    resamp_grid = [1, 5, 9, 13, 14, 17, 21, 22, 23] * u.AA

    inst = SplineInterpolatedResampler()
    results = inst(input_spectra, resamp_grid)

    assert_quantity_allclose(results.flux,
                             np.array([np.nan, 3.98808594, 6.94042969, 6.45869141,
                                       5.89921875, 7.29736328, np.nan, np.nan,
                                       np.nan])*u.mJy)

    inst = SplineInterpolatedResampler(extrapolation_treatment="truncate")
    input_spectra = Spectrum1D(spectral_axis=wave_val * u.AA, flux=flux_val * u.mJy)
    results = inst(input_spectra, resamp_grid)

    assert_quantity_allclose(results.flux,
                             np.array([3.98808594, 6.94042969, 6.45869141,
                                       5.89921875, 7.29736328])*u.mJy)


@pytest.mark.parametrize("edgetype,lastvalue",
                         [("nan_fill", np.nan), ("zero_fill", 0)])
def test_resample_edges(edgetype, lastvalue, all_resamplers):
    input_spectrum = Spectrum1D(spectral_axis=[2, 4, 12, 16, 20] * u.micron,
                                flux=[1, 3, 7, 6, 20] * u.mJy)
    resamp_grid = [1, 3, 7, 16, 20, 100] * u.micron

    resampler = all_resamplers(edgetype)
    resampled = resampler(input_spectrum, resamp_grid)
    if lastvalue is np.nan:
        assert np.isnan(resampled.flux[-1])
    else:
        assert resampled.flux[-1] == lastvalue


def test_resample_different_units(all_resamplers):
    input_spectrum = Spectrum1D(spectral_axis=[5000, 6000, 7000] * u.AA,
                                flux=[1, 2, 3] * u.mJy)
    resampler = all_resamplers("nan_fill")
    if all_resamplers == FluxConservingResampler:
        pytest.xfail('flux conserving resampler cannot yet handle differing units')

    resamp_grid = [5500, 6500]*u.nm
    resampled = resampler(input_spectrum, resamp_grid)
    assert np.all(np.isnan(resampled.flux))

    resamp_grid = [550, 650]*u.nm
    resampled = resampler(input_spectrum, resamp_grid)
    assert not np.any(np.isnan(resampled.flux))

    resamp_grid = [550, 650]*u.nm
    resampled = resampler(input_spectrum, resamp_grid)

    # Test conversion to velocity grid
    rest_wavelength = 656.2 * u.nm
    wavelengths = np.linspace(640, 672, 10) * u.nm
    flux = np.ones(10) * u.mJy
    spec1d = Spectrum1D(spectral_axis=wavelengths, velocity_convention="optical", flux=flux)
    spec1d.spectral_axis.doppler_rest = rest_wavelength

    velocities = np.linspace(-1000, 1000, 5) * u.km/u.s
    velocity_grid = SpectralAxis(velocities, doppler_rest=rest_wavelength,
                                 doppler_convention="optical")
    velocity_binned = resampler(spec1d, velocity_grid)
    assert not np.any(np.isnan(velocity_binned.flux))


def test_resample_uncs(all_resamplers):
    sdunc = StdDevUncertainty([0.1, 0.2, 0.3]*u.mJy)
    input_spectrum = Spectrum1D(spectral_axis=[5000, 6000, 7000] * u.AA,
                                flux=[1, 2, 3] * u.mJy,
                                uncertainty=sdunc)

    resampled = all_resamplers()(input_spectrum, [5500, 6500]*u.AA)
    if all_resamplers == FluxConservingResampler:
        # special-cased because it switches the unc to inverse variance by construction
        assert resampled.uncertainty.unit == sdunc.unit**-2
    else:
        assert resampled.uncertainty.unit == sdunc.unit
        assert resampled.uncertainty.uncertainty_type == sdunc.uncertainty_type


def test_fluxconservingresampler_against_spectres():
    # this test compares the output from FluxConservingResampler to the output
    # from running the code in https://github.com/ACCarnall/SpectRes, which
    # contains the algorithm that FluxConservedResampler is also
    # implementing (slightly differently, but results should always agree)

    spectres_fluxes = np.array([2., 3., 3.66666667, 6.3125])
    spectres_errs = np.array([0.25, 0.15, 0.30413813, 0.49193694])

    flux = [1, 2, 3, 4, 5, 6, 7]*u.AB
    spectral_axis = [1, 2, 3, 7, 9, 10, 14]*u.AA
    new_spec_axis = [2, 3, 4, 12]*u.AA
    errs = VarianceUncertainty([0.5, 0.25, 0.15, 0.45, 0.75, 1.5, 0.1])
    input_spectra = Spectrum1D(spectral_axis=spectral_axis, flux=flux,
                               uncertainty=errs)

    inst = FluxConservingResampler()
    fluxc_resampled = inst(input_spectra, new_spec_axis)
    fluxc_flux = fluxc_resampled.flux.value
    # have to do 1 / errs because InverseVariance is always returned
    fluxc_errs = 1 / fluxc_resampled.uncertainty.array

    assert_quantity_allclose(fluxc_flux, spectres_fluxes)
    assert_quantity_allclose(fluxc_errs, spectres_errs)
