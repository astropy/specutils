import astropy.units as u
import numpy as np
import pytest

from ..spectra.spectrum1d import Spectrum1D
from .spectral_examples import simulated_spectra


def test_spectral_axes():
    flux1 = (np.random.sample(49) * 100).astype(int)
    flux2 = (np.random.sample(49) * 100).astype(int)

    flux3 = flux1 + flux2

    spec1 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                       flux=flux1 * u.Jy)
    spec2 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                       flux=flux2 * u.Jy)

    spec3 = spec1 + spec2

    assert np.allclose(spec3.flux.value, flux3)


def test_add_basic_spectra(simulated_spectra):

    # Get the numpy array of data
    flux1 = simulated_spectra.s1_um_mJy_e1_flux
    flux2 = simulated_spectra.s1_um_mJy_e2_flux
    flux3 = flux1 + flux2

    # Calculate using the spectrum1d/nddata code
    spec3 = simulated_spectra.s1_um_mJy_e1 + simulated_spectra.s1_um_mJy_e2

    assert np.allclose(spec3.flux.value, flux3)


def test_add_diff_flux_prefix(simulated_spectra):

    # Get the numpy array of data
    # this assumes output will be in mJy units
    flux1 = simulated_spectra.s1_AA_mJy_e3_flux
    flux2 = simulated_spectra.s1_AA_nJy_e4_flux
    flux3 = flux1 + (flux2 / 1000000)

    # Calculate using the spectrum1d/nddata code
    spec3 = simulated_spectra.s1_AA_mJy_e3 + simulated_spectra.s1_AA_nJy_e4

    assert np.allclose(spec3.flux.value, flux3)


def test_subtract_basic_spectra(simulated_spectra):

    # Get the numpy array of data
    flux1 = simulated_spectra.s1_um_mJy_e1_flux
    flux2 = simulated_spectra.s1_um_mJy_e2_flux
    flux3 = flux2 - flux1

    # Calculate using the spectrum1d/nddata code
    spec3 = simulated_spectra.s1_um_mJy_e2 - simulated_spectra.s1_um_mJy_e1

    assert np.allclose(spec3.flux.value, flux3)


def test_divide_basic_spectra(simulated_spectra):

    # Get the numpy array of data
    flux1 = simulated_spectra.s1_um_mJy_e1_flux
    flux2 = simulated_spectra.s1_um_mJy_e2_flux
    flux3 = flux1 / flux2

    # Calculate using the spectrum1d/nddata code
    spec3 = simulated_spectra.s1_um_mJy_e1 / simulated_spectra.s1_um_mJy_e2

    assert np.allclose(spec3.flux.value, flux3)


def test_multiplication_basic_spectra(simulated_spectra):

    # Get the numpy array of data
    flux1 = simulated_spectra.s1_um_mJy_e1_flux
    flux2 = simulated_spectra.s1_um_mJy_e2_flux
    flux3 = flux1 * flux2

    # Calculate using the spectrum1d/nddata code
    spec3 = simulated_spectra.s1_um_mJy_e1 * simulated_spectra.s1_um_mJy_e2

    assert np.allclose(spec3.flux.value, flux3)


@pytest.mark.xfail(raises=ValueError)
def test_add_diff_spectral_axis(simulated_spectra):

    # Calculate using the spectrum1d/nddata code
    spec3 = simulated_spectra.s1_um_mJy_e1 + simulated_spectra.s1_AA_mJy_e3
