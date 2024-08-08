import astropy.units as u
from astropy.tests.helper import assert_quantity_allclose
import numpy as np

from ..spectra.spectrum1d import Spectrum1D


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
    spec3 = simulated_spectra.s1_um_mJy_e1 + simulated_spectra.s1_AA_mJy_e3  # noqa


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

    assert spec3.mask[nan_idx].all() == True  # noqa


def test_with_constants(simulated_spectra):
    spec = simulated_spectra.s1_um_mJy_e1

    # Test that doing arithmetic with a constant to the right of the spectrum succeeds
    assert_quantity_allclose((2 * spec).flux, spec.flux * 2)

    r_add_result = 2 + spec
    l_add_result = spec + 2
    assert_quantity_allclose(r_add_result.flux, l_add_result.flux)

    r_sub_result = 2 - spec
    l_sub_result = -1 * (spec - 2)
    assert_quantity_allclose(r_sub_result.flux, l_sub_result.flux)


def test_arithmetic_after_shift(simulated_spectra):
    spec = simulated_spectra.s1_um_mJy_e1
    spec.shift_spectrum_to(redshift = 1)

    # Test that doing arithmetic preserves the shifted spectral axis
    spec *= 2
    assert_quantity_allclose(spec.spectral_axis, 2*np.linspace(0.4, 1.05, 100)*u.um)
