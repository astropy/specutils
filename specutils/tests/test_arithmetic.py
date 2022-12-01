from copy import deepcopy

import astropy.units as u
import numpy as np
import pytest
from astropy.tests.helper import assert_quantity_allclose

from specutils.spectra.spectrum1d import Spectrum1D


# Mimic specreduce background object a little.
class _MockBackground:
    def __init__(self, spec):
        self.image = spec

    def __rsub__(self, other):
        return 42  # Does not matter what, just want to make sure this is called.


class TestMathWithAllOnes:
    def setup_class(self):
        flux = np.ones(10) * u.nJy
        wave = (np.arange(flux.size) + 1) * u.um
        self.spec = Spectrum1D(spectral_axis=wave, flux=flux)

    @pytest.mark.parametrize('wave_unit', ('same', 'Angstrom'))
    def test_math_with_spectral_axes_same(self, wave_unit):
        if wave_unit == 'same':
            new_spec = self.spec
        else:  # Same spectral axis but in different unit.
            new_spec = Spectrum1D(spectral_axis=self.spec.spectral_axis.to(u.Unit(wave_unit)),
                                  flux=self.spec.flux)

        spec_added = self.spec + new_spec
        assert_quantity_allclose(spec_added.flux, 2 * u.nJy)

        spec_subbed = self.spec - new_spec
        assert_quantity_allclose(spec_subbed.flux, 0 * u.nJy)

        spec_mul = self.spec * self.spec
        assert_quantity_allclose(spec_mul.flux, 1 * (u.nJy * u.nJy))

        spec_div = self.spec / self.spec
        assert_quantity_allclose(spec_div.flux, 1)

    def test_math_with_spectral_axes_different(self):
        new_wave = self.spec.spectral_axis + (1 * u.um)
        new_spec = Spectrum1D(spectral_axis=new_wave, flux=self.spec.flux)

        with pytest.raises(ValueError, match='Mismatched spectral_axis'):
            self.spec + new_spec

        with pytest.raises(ValueError, match='Mismatched spectral_axis'):
            self.spec - new_spec

        with pytest.raises(ValueError, match='Mismatched spectral_axis'):
            self.spec * new_spec

        with pytest.raises(ValueError, match='Mismatched spectral_axis'):
            self.spec / new_spec

    def test_math_with_plain_number(self):
        spec_added_fwd = self.spec + 1
        assert_quantity_allclose(spec_added_fwd.flux, 2 * u.nJy)

        spec_added_bak = 1 + self.spec
        assert_quantity_allclose(spec_added_bak.flux, 2 * u.nJy)

        spec_subbed_fwd = self.spec - 1
        assert_quantity_allclose(spec_subbed_fwd.flux, 0 * u.nJy)

        spec_subbed_bak = 2 - self.spec
        assert_quantity_allclose(spec_subbed_bak.flux, 1 * u.nJy)

        spec_mul_fwd = self.spec * 2
        assert_quantity_allclose(spec_mul_fwd.flux, 2 * u.nJy)

        spec_mul_bak = 2 * self.spec
        assert_quantity_allclose(spec_mul_bak.flux, 2 * u.nJy)

        spec_div = self.spec / 2
        assert_quantity_allclose(spec_div.flux, 0.5 * u.nJy)

        with pytest.raises(TypeError, match='unsupported operand'):
            1 / self.spec

        with pytest.raises(NotImplementedError, match='Cannot operate on ndarray class'):
            self.spec * np.ones(self.spec.flux.shape)

    def test_math_with_quantity(self):
        q = 1 * u.uJy  # 1000 nJy

        spec_added_fwd = self.spec + q
        assert_quantity_allclose(spec_added_fwd.flux, 1001 * u.nJy)

        # Quantity has its own __add__
        with pytest.raises(ValueError, match='Value not scalar compatible'):
            q + self.spec

        spec_subbed_fwd = self.spec - q
        assert_quantity_allclose(spec_subbed_fwd.flux, -999 * u.nJy)

        # Quantity has its own __sub__
        with pytest.raises(ValueError, match='Value not scalar compatible'):
            q - self.spec

        spec_mul_fwd = self.spec * q
        assert_quantity_allclose(spec_mul_fwd.flux, 1000 * (u.nJy * u.nJy))

        # Quantity has its own __mul__
        spec_mul_bak = q * self.spec
        assert isinstance(spec_mul_bak, u.Quantity)  # Becomes Quantity

        spec_div = self.spec / q
        assert_quantity_allclose(spec_div.flux, 0.001)

        # Quantity has its own __truediv__
        with pytest.raises(TypeError, match='unsupported operand'):
            q / self.spec

        with pytest.raises(ValueError, match='Quantity must be scalar'):
            self.spec * (self.spec.flux)

    def test_math_with_ndcube(self):
        from astropy.wcs import WCS
        from ndcube import NDCube

        data = np.ones((4, 4, 10)) * u.nJy
        wcs = WCS(naxis=3)
        wcs.wcs.ctype = 'WAVE', 'RA--TAN', 'DEC-TAN'
        wcs.wcs.cunit = 'Angstrom', 'deg', 'deg'
        wcs.wcs.cdelt = 0.2, 0.5, 0.4
        wcs.wcs.crpix = 0, 2, 2
        wcs.wcs.crval = 10, 0.5, 1
        wcs.wcs.cname = 'wavelength', 'lon', 'lat'

        ndc = NDCube(data, wcs=wcs)
        spec3d = Spectrum1D(flux=data, wcs=wcs)

        spec_added = spec3d + ndc
        assert_quantity_allclose(spec_added.flux, 2 * u.nJy)

        spec_subbed = spec3d - ndc
        assert_quantity_allclose(spec_subbed.flux, 0 * u.nJy)

        spec_mul = spec3d * ndc
        assert_quantity_allclose(spec_mul.flux, 1 * (u.nJy * u.nJy))

        spec_div = spec3d / ndc
        assert_quantity_allclose(spec_div.flux, 1)

        # Also test 3D vs 1D Spectrum1D
        with pytest.raises(ValueError, match='Mismatched spectral_axis'):
            spec3d + self.spec

    def test_nd_vs_1d(self):
        spec2d = Spectrum1D(flux=np.ones((10, 10)) * u.nJy, spectral_axis=self.spec.spectral_axis)
        spec_added_2d = spec2d + self.spec
        assert_quantity_allclose(spec_added_2d.flux, 2 * u.nJy)

        spec3d = Spectrum1D(flux=np.ones((10, 10, 10)) * u.nJy, spectral_axis=self.spec.spectral_axis)
        spec_added_3d = spec3d + self.spec
        assert_quantity_allclose(spec_added_3d.flux, 2 * u.nJy)

    def test_mask_nans(self):
        new_flux = deepcopy(self.spec.flux)
        nan_idx = [1, 3, 5]
        new_flux[nan_idx] = np.nan
        new_spec = Spectrum1D(spectral_axis=self.spec.spectral_axis, flux=new_flux)
        spec_added = self.spec + new_spec
        assert spec_added.mask[nan_idx].all()

    def test_specreduce_bg_rsub(self):
        bg = _MockBackground(self.spec)
        assert (self.spec - bg) == 42


def test_add_basic_spectra(simulated_spectra):

    # Get the numpy array of data
    flux1 = simulated_spectra.s1_um_mJy_e1_flux
    flux2 = simulated_spectra.s1_um_mJy_e2_flux
    flux3 = flux1 + flux2

    # Calculate using the spectrum1d/nddata code
    spec3 = simulated_spectra.s1_um_mJy_e1 + simulated_spectra.s1_um_mJy_e2

    np.testing.assert_allclose(spec3.flux.value, flux3)


def test_add_diff_flux_prefix(simulated_spectra):

    # Get the numpy array of data
    # this assumes output will be in mJy units
    flux1 = simulated_spectra.s1_AA_mJy_e3_flux
    flux2 = simulated_spectra.s1_AA_nJy_e4_flux
    flux3 = flux1 + (flux2 / 1000000)

    # Calculate using the spectrum1d/nddata code
    spec3 = simulated_spectra.s1_AA_mJy_e3 + simulated_spectra.s1_AA_nJy_e4

    np.testing.assert_allclose(spec3.flux.value, flux3)


def test_subtract_basic_spectra(simulated_spectra):

    # Get the numpy array of data
    flux1 = simulated_spectra.s1_um_mJy_e1_flux
    flux2 = simulated_spectra.s1_um_mJy_e2_flux
    flux3 = flux2 - flux1

    # Calculate using the spectrum1d/nddata code
    spec3 = simulated_spectra.s1_um_mJy_e2 - simulated_spectra.s1_um_mJy_e1

    np.testing.assert_allclose(spec3.flux.value, flux3)


def test_divide_basic_spectra(simulated_spectra):

    # Get the numpy array of data
    flux1 = simulated_spectra.s1_um_mJy_e1_flux
    flux2 = simulated_spectra.s1_um_mJy_e2_flux
    flux3 = flux1 / flux2

    # Calculate using the spectrum1d/nddata code
    spec3 = simulated_spectra.s1_um_mJy_e1 / simulated_spectra.s1_um_mJy_e2

    np.testing.assert_allclose(spec3.flux.value, flux3)


def test_multiplication_basic_spectra(simulated_spectra):

    # Get the numpy array of data
    flux1 = simulated_spectra.s1_um_mJy_e1_flux
    flux2 = simulated_spectra.s1_um_mJy_e2_flux
    flux3 = flux1 * flux2

    # Calculate using the spectrum1d/nddata code
    spec3 = simulated_spectra.s1_um_mJy_e1 * simulated_spectra.s1_um_mJy_e2

    np.testing.assert_allclose(spec3.flux.value, flux3)


def test_masks(simulated_spectra):
    masked_spec = simulated_spectra.s1_um_mJy_e1_masked

    masked_sum = masked_spec + masked_spec
    assert np.all(masked_sum.mask == simulated_spectra.s1_um_mJy_e1_masked.mask)

    masked_sum.mask[:50] = True
    masked_diff = masked_sum - masked_spec
    assert_quantity_allclose(masked_diff.flux, masked_spec.flux)
    assert np.all(masked_diff.mask == masked_sum.mask | masked_spec.mask)


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
