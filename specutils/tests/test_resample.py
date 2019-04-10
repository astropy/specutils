import numpy as np
import pytest
import astropy.units as u

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
