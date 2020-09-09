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


def test_add_diff_spectral_axis(simulated_spectra):

    # Calculate using the spectrum1d/nddata code
    spec3 = simulated_spectra.s1_um_mJy_e1 + simulated_spectra.s1_AA_mJy_e3


def test_masks(simulated_spectra):
    masked_spec = simulated_spectra.s1_um_mJy_e1_masked

    masked_sum = masked_spec + masked_spec
    assert np.all(masked_sum.mask == simulated_spectra.s1_um_mJy_e1_masked.mask)

    masked_sum.mask[:50] = True
    masked_diff = masked_sum - masked_spec
    assert u.allclose(masked_diff.flux, masked_spec.flux)
    assert np.all(masked_diff.mask == masked_sum.mask | masked_spec.mask)


def test_mask_nans():
    flux1 = np.random.random(10)
    flux2 = np.random.random(10)
    nan_idx = [1, 3, 5]
    flux2[nan_idx] = np.nan
    spec1 = Spectrum1D(spectral_axis=np.arange(10) * u.nm, flux=flux1 * u.Jy)
    spec2 = Spectrum1D(spectral_axis=np.arange(10) * u.nm, flux=flux2 * u.Jy)

    spec3 = spec1 + spec2

    assert spec3.mask[nan_idx].all() == True
