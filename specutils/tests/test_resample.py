import numpy as np
import pytest
import astropy.units as u
from astropy.tests.helper import assert_quantity_allclose

from ..spectra.spectrum1d import Spectrum1D
from ..tests.spectral_examples import simulated_spectra
from ..manipulation.resample import FluxConservingResample

# todo: Should add tests for different weighting options once those
# are more solidified.


def test_same_grid_fluxconserving(simulated_spectra):
    """
    Test that feeding in the original dispersion axis returns the
    same flux after resampling.
    """

    inst = FluxConservingResample()
    results = inst(simulated_spectra.s1_um_mJy_e1,
                   simulated_spectra.s1_um_mJy_e1.wavelength, weights=None)

    assert  np.allclose(np.array(simulated_spectra.s1_um_mJy_e1.flux),
                        np.array(results.flux))


def test_expanded_grid_fluxconserving():
    """
    New dispersion axis has more bins then input dispersion axis
    """

    flux_val = np.array([1, 3, 7, 6, 20])
    wave_val = np.array([2, 4, 12, 16, 20])
    input_spectra = Spectrum1D(flux_val*u.mJy, wave_val* u.AA)
    resamp_grid = np.array([1, 5, 9, 13, 14, 17, 21, 22, 23])

    inst = FluxConservingResample()
    results = inst(input_spectra, resamp_grid, weights=None)

    assert results.flux.unit == u.mJy
    assert np.allclose(np.array(results.flux),
                        np.array([0., 3., 6.13043478, 7., 6.33333333, 10.,
                                  20., 0., 0.]))


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


@pytest.mark.parametrize("specflux,specwave,outwave", [
    ([1, 3, 2], [4000, 5000, 6000], np.linspace(4000, 6000, 4)),
    ([1, 3, 2,1], np.linspace(4000, 6000, 4), [4000, 5000, 6000])
    ])
def test_flux_conservation(specflux, specwave, outwave):
    """
    A few simple cases to programatically ensure flux is conserved in the
    resampling algorithm
    """
    in_spec = Spectrum1D(spectral_axis=specwave*u.AA, flux=specflux*u.AB)
    out_spec = FluxConservingResample()(in_spec, outwave*u.AA, weights=None)

    in_dwl = delta_wl(in_spec.spectral_axis)
    out_dwl = delta_wl(out_spec.spectral_axis)

    assert assert_quantity_allclose(np.sum(in_spec.flux * in_dwl),
                                    np.sum(out_spec.flux * out_dwl))
