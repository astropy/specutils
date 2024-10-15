import numpy as np
import pytest

import astropy.units as u
from astropy.nddata import StdDevUncertainty
from astropy.tests.helper import assert_quantity_allclose

from ..spectra import Spectrum1D, SpectralRegion
from ..manipulation import extract_region, extract_bounding_spectral_region, spectral_slab
from ..manipulation.utils import linear_exciser

FLUX_ARRAY = [1605.71612173, 1651.41650744, 2057.65798618, 2066.73502361, 1955.75832537,
              1670.52711471, 1491.10034446, 1637.08084112, 1471.28982259, 1299.19484483,
              1423.11195734, 1226.74494917, 1572.31888312, 1311.50503403, 1474.05051673,
              1335.39944397, 1420.61880528, 1433.18623759, 1290.26966668, 1605.67341284,
              1528.52281708, 1592.74392861, 1568.74162534, 1435.29407808, 1536.68040935,
              1157.33825995, 1136.12679394, 999.92394692, 1038.61546167, 1011.60297294]


def test_region_simple(simulated_spectra):

    spectrum = simulated_spectra.s1_um_mJy_e1
    uncertainty = StdDevUncertainty(np.full(spectrum.flux.shape, 0.1) * u.mJy)
    spectrum.uncertainty = uncertainty

    region = SpectralRegion(0.6*u.um, 0.8*u.um)

    sub_spectrum = extract_region(spectrum, region)

    sub_spectrum_flux_expected = np.array(FLUX_ARRAY)

    assert_quantity_allclose(sub_spectrum.flux.value, sub_spectrum_flux_expected)


def test_pixel_spectralaxis_extraction():
    """
    Tests a region extraction on a Spectrum1D with a u.pix defined spectral axis
    """
    flux_unit = u.dimensionless_unscaled
    spec_unit = u.pix

    spec1d = Spectrum1D(spectral_axis=np.arange(5100, 5300)*spec_unit,
                        flux=np.random.randn(200)*flux_unit)

    # Case 1: Region is safely within the bounds of a continuously-defined Spec1D
    # Region is intentionally defined "between values" to test for index rounding
    region = SpectralRegion.from_center(center=5200.5*spec_unit,
                                        width=100*spec_unit)
    extracted_spec1d = extract_region(spec1d, region)

    # Extraction is of the right shape and value
    assert extracted_spec1d.shape == (101,)
    assert_quantity_allclose(extracted_spec1d.spectral_axis,
                             spec1d.spectral_axis[50:151])
    assert_quantity_allclose(extracted_spec1d.flux,
                             spec1d.flux[50:151])

    # Follow Python slicing conventions:
    # Lower slice is inclusive
    assert extracted_spec1d.spectral_axis[0].quantity <= region.lower
    # Upper slice is exclusive
    assert extracted_spec1d.spectral_axis[-1].quantity < region.upper

    # Case 2: Region is outside the lower bounds of the Spec1D
    region2 = SpectralRegion.from_center(center=spec1d.spectral_axis[0].quantity,
                                         width=100*spec_unit)
    extracted_spec1d_2 = extract_region(spec1d, region2)

    assert extracted_spec1d_2.shape == (50,)  # Included upper bound is exclusive
    assert_quantity_allclose(extracted_spec1d_2.spectral_axis,
                             spec1d.spectral_axis[0:50])
    assert_quantity_allclose(extracted_spec1d_2.flux,
                             spec1d.flux[0:50])

    # Case 3: Region is outside the upper bounds of the Spec1D
    region3 = SpectralRegion.from_center(center=spec1d.spectral_axis[-1].quantity,
                                         width=100*spec_unit)
    extracted_spec1d_3 = extract_region(spec1d, region3)

    assert extracted_spec1d_3.shape == (51,)  # Included lower bound is inclusive (+1)
    assert_quantity_allclose(extracted_spec1d_3.spectral_axis,
                             spec1d.spectral_axis[149:])
    assert_quantity_allclose(extracted_spec1d_3.flux,
                             spec1d.flux[149:])

    # Case 4: Compound region with the two definitions above
    extracted_spec1d_4 = extract_region(spec1d, (region2+region3))

    # Each extracted part of compound region should be identical to the case 2 and 3 extractions
    # If they're identical, then the previous value checks for case 2/3 should cover this as well
    assert len(extracted_spec1d_4) == 2

    for original_case_spectra, compound_spectra in [(extracted_spec1d_2, extracted_spec1d_4[0]),
                                                    (extracted_spec1d_3, extracted_spec1d_4[1])]:
        assert original_case_spectra.shape == compound_spectra.shape
        assert_quantity_allclose(original_case_spectra.spectral_axis,
                                 compound_spectra.spectral_axis)
        assert_quantity_allclose(original_case_spectra.flux,
                                 compound_spectra.flux)

    # Case 5: Region is entirely outside bounds of the spectra should return nothing
    upper_bound = spec1d.spectral_axis[-1].quantity
    region = SpectralRegion((upper_bound + 10*spec_unit), (upper_bound + 100*spec_unit))
    extracted_spec1d = extract_region(spec1d, region)
    assert extracted_spec1d.shape == (0,)

    lower_bound = spec1d.spectral_axis[0].quantity
    region = SpectralRegion((lower_bound - 100*spec_unit), (lower_bound - 10*spec_unit))
    extracted_spec1d = extract_region(spec1d, region)
    assert extracted_spec1d.shape == (0,)


def test_slab_simple(simulated_spectra):
    spectrum = simulated_spectra.s1_um_mJy_e1
    uncertainty = StdDevUncertainty(np.full(spectrum.flux.shape, 0.1) * u.mJy)
    spectrum.uncertainty = uncertainty

    sub_spectrum = spectral_slab(spectrum, 0.6*u.um, 0.8*u.um)

    sub_spectrum_flux_expected = np.array(FLUX_ARRAY)

    assert_quantity_allclose(sub_spectrum.flux.value, sub_spectrum_flux_expected)


def test_slab_pixels():
    range_iter = range(5)
    spectrum = Spectrum1D(flux=np.stack(
        [np.zeros((2, 2)) + i for i in range_iter], axis=-1) * u.nJy)
    for i in range_iter:
        sub_spectrum = spectral_slab(spectrum, i * u.pix, (i + 0.5) * u.pix)
        assert_quantity_allclose(sub_spectrum.flux, i * u.nJy)
        assert_quantity_allclose(sub_spectrum.spectral_axis, i * u.pix)


def test_region_ghz(simulated_spectra):
    spectrum = Spectrum1D(flux=simulated_spectra.s1_um_mJy_e1.flux,
                          spectral_axis=simulated_spectra.s1_um_mJy_e1.frequency)

    region = SpectralRegion(499654.09666667*u.GHz, 374740.5725*u.GHz)

    sub_spectrum = extract_region(spectrum, region)

    sub_spectrum_flux_expected = FLUX_ARRAY * u.mJy

    assert_quantity_allclose(sub_spectrum.flux, sub_spectrum_flux_expected)


def test_region_ghz_spectrum_wave(simulated_spectra):
    spectrum = simulated_spectra.s1_um_mJy_e1
    region = SpectralRegion(499654.09666667*u.GHz, 374740.5725*u.GHz)

    sub_spectrum = extract_region(spectrum, region)

    sub_spectrum_flux_expected = FLUX_ARRAY * u.mJy
    assert_quantity_allclose(sub_spectrum.flux, sub_spectrum_flux_expected)


def test_region_simple_check_ends(simulated_spectra):
    np.random.seed(42)

    spectrum = Spectrum1D(spectral_axis=np.linspace(1, 25, 25)*u.um, flux=np.random.random(25)*u.Jy)
    region = SpectralRegion(8*u.um, 15*u.um)

    sub_spectrum = extract_region(spectrum, region)

    assert sub_spectrum.spectral_axis.value[0] == 8
    assert sub_spectrum.spectral_axis.value[-1] == 15

    region = SpectralRegion(0*u.um, 15*u.um)
    sub_spectrum = extract_region(spectrum, region)
    assert sub_spectrum.spectral_axis.value[0] == 1

    region = SpectralRegion(8*u.um, 30*u.um)
    sub_spectrum = extract_region(spectrum, region)
    assert sub_spectrum.spectral_axis.value[-1] == 25


def test_region_empty():
    empty_spectrum = Spectrum1D(spectral_axis=[]*u.um, flux=[]*u.Jy)

    # Region past upper range of spectrum
    spectrum = Spectrum1D(spectral_axis=np.linspace(1, 25, 25)*u.um, flux=np.ones(25)*u.Jy)
    region = SpectralRegion(28*u.um, 30*u.um)
    sub_spectrum = extract_region(spectrum, region)

    assert np.allclose(sub_spectrum.spectral_axis.value, empty_spectrum.spectral_axis.value)
    assert sub_spectrum.spectral_axis.unit == empty_spectrum.spectral_axis.unit

    assert np.allclose(sub_spectrum.flux.value, empty_spectrum.flux.value)
    assert sub_spectrum.flux.unit == empty_spectrum.flux.unit

    # Region below lower range of spectrum
    spectrum = Spectrum1D(spectral_axis=np.linspace(1, 25, 25)*u.um, flux=np.ones(25)*u.Jy)
    region = SpectralRegion(0.1*u.um, 0.3*u.um)
    sub_spectrum = extract_region(spectrum, region)

    assert np.allclose(sub_spectrum.spectral_axis.value, empty_spectrum.spectral_axis.value)
    assert sub_spectrum.spectral_axis.unit == empty_spectrum.spectral_axis.unit

    assert np.allclose(sub_spectrum.flux.value, empty_spectrum.flux.value)
    assert sub_spectrum.flux.unit == empty_spectrum.flux.unit

    # Region below lower range of spectrum and upper range in the spectrum.
    spectrum = Spectrum1D(spectral_axis=np.linspace(1, 25, 25)*u.um, flux=2*np.linspace(1, 25, 25)*u.Jy)
    region = SpectralRegion(0.1*u.um, 3.3*u.um)
    sub_spectrum = extract_region(spectrum, region)

    assert np.allclose(sub_spectrum.spectral_axis.value, [1, 2, 3])
    assert sub_spectrum.spectral_axis.unit == empty_spectrum.spectral_axis.unit

    assert np.allclose(sub_spectrum.flux.value, [2, 4, 6])
    assert sub_spectrum.flux.unit == empty_spectrum.flux.unit

    # Region has lower and upper bound the same
    with pytest.raises(Exception):
        region = SpectralRegion(3*u.um, 3*u.um)


def test_region_descending(simulated_spectra):
    spectrum = simulated_spectra.s1_um_mJy_e1
    uncertainty = StdDevUncertainty(np.full(spectrum.flux.shape, 0.1) * u.mJy)
    spectrum.uncertainty = uncertainty

    region = SpectralRegion(0.8*u.um, 0.6*u.um)

    sub_spectrum = extract_region(spectrum, region)

    sub_spectrum_flux_expected = np.array(FLUX_ARRAY)

    assert_quantity_allclose(sub_spectrum.flux.value, sub_spectrum_flux_expected)


def test_descending_spectral_axis(simulated_spectra):
    spectrum = simulated_spectra.s1_um_mJy_e1_desc

    sub_spectrum_flux_expected = np.array(FLUX_ARRAY[::-1])

    region = SpectralRegion(0.8*u.um, 0.6*u.um)
    sub_spectrum = extract_region(spectrum, region)

    assert_quantity_allclose(sub_spectrum.flux.value, sub_spectrum_flux_expected)

    region = SpectralRegion(0.6*u.um, 0.8*u.um)
    sub_spectrum = extract_region(spectrum, region)

    assert_quantity_allclose(sub_spectrum.flux.value, sub_spectrum_flux_expected)


def test_region_two_sub(simulated_spectra):
    spectrum = simulated_spectra.s1_um_mJy_e1
    uncertainty = StdDevUncertainty(np.full(spectrum.flux.shape, 0.1) * u.mJy)
    spectrum.uncertainty = uncertainty

    region = SpectralRegion([(0.6*u.um, 0.8*u.um), (0.86*u.um, 0.89*u.um)])

    sub_spectra = extract_region(spectrum, region)

    # Confirm the end points of the subspectra are correct
    assert_quantity_allclose(sub_spectra[0].spectral_axis[[0, -1]],
                       [0.6035353535353536, 0.793939393939394]*u.um)

    assert_quantity_allclose(sub_spectra[1].spectral_axis[[0, -1]],
                       [0.8661616161616162, 0.8858585858585859]*u.um)

    sub_spectrum_0_flux_expected = FLUX_ARRAY * u.mJy

    sub_spectrum_1_flux_expected = [1337.65312465, 1263.48914109,
                                    1589.81797876, 1548.46068415]*u.mJy

    assert_quantity_allclose(sub_spectra[0].flux, sub_spectrum_0_flux_expected)
    assert_quantity_allclose(sub_spectra[1].flux, sub_spectrum_1_flux_expected)

    # also ensure this works if the multi-region is expressed as a single
    # Quantity
    region2 = SpectralRegion([(0.6, 0.8), (0.86, 0.89)]*u.um)
    sub_spectra2 = extract_region(spectrum, region2)
    assert_quantity_allclose(sub_spectra[0].flux, sub_spectra2[0].flux)
    assert_quantity_allclose(sub_spectra[1].flux, sub_spectra2[1].flux)

    # Check that the return_single_spectrum argument works properly
    concatenated_spectrum = extract_region(spectrum, region2, return_single_spectrum=True)
    assert concatenated_spectrum.flux.shape == (34,)
    assert np.all(concatenated_spectrum.flux[0:30] == sub_spectra2[0].flux)
    assert np.all(concatenated_spectrum.flux[30:34] == sub_spectra2[1].flux)


def test_bounding_region(simulated_spectra):
    spectrum = simulated_spectra.s1_um_mJy_e1
    uncertainty = StdDevUncertainty(np.full(spectrum.flux.shape, 0.1) * u.mJy)
    spectrum.uncertainty = uncertainty

    region = SpectralRegion([(0.6*u.um, 0.8*u.um), (0.86*u.um, 0.89*u.um)])

    extracted_spectrum = extract_bounding_spectral_region(spectrum, region)

    # Confirm the end points are correct
    assert_quantity_allclose(extracted_spectrum.spectral_axis[[0, -1]],
                       [0.6035353535353536, 0.8858585858585859]*u.um)

    flux_expected = (FLUX_ARRAY + [948.81864554, 1197.84859443, 1069.75268943,
                                   1118.27269184, 1301.7695563, 1206.62880648,
                                   1518.16549319, 1256.84259015, 1638.76791267,
                                   1562.05642302, 1337.65312465, 1263.48914109,
                                   1589.81797876, 1548.46068415])*u.mJy

    assert_quantity_allclose(extracted_spectrum.flux, flux_expected)

    # also ensure this works if the multi-region is expressed as a single
    # Quantity
    region2 = SpectralRegion([(0.6, 0.8), (0.86, 0.89)]*u.um)
    extracted_spectrum2 = extract_bounding_spectral_region(spectrum, region2)
    assert_quantity_allclose(extracted_spectrum2.spectral_axis[[0, -1]],
                       [0.6035353535353536, 0.8858585858585859]*u.um)
    assert_quantity_allclose(extracted_spectrum2.flux, flux_expected)


def test_extract_region_pixels():
    spectrum = Spectrum1D(spectral_axis=np.linspace(4000, 10000, 25)*u.AA,
                          flux=np.arange(25)*u.Jy)
    region = SpectralRegion(10*u.pixel, 12*u.pixel)

    extracted = extract_region(spectrum, region)

    assert_quantity_allclose(extracted.flux, [10, 11]*u.Jy)


def test_extract_region_mismatched_units():
    spectrum = Spectrum1D(spectral_axis=np.arange(25)*u.nm,
                          flux=np.arange(25)*u.Jy)

    region = SpectralRegion(100*u.AA, 119*u.AA)

    extracted = extract_region(spectrum, region)

    assert_quantity_allclose(extracted.flux, [10, 11]*u.Jy)


@pytest.mark.filterwarnings('ignore:A SpectralRegion with multiple subregions was provided')
def test_linear_excise_invert_from_spectrum():
    spec = Spectrum1D(flux=np.random.sample(100) * u.Jy,
                      spectral_axis=np.arange(100) * u.AA)
    inc_regs = (SpectralRegion(0 * u.AA, 50 * u.AA) +
                SpectralRegion(60 * u.AA, 80 * u.AA) +
                SpectralRegion(90 * u.AA, 110 * u.AA))
    exc_regs = inc_regs.invert_from_spectrum(spec)

    excised_spec = linear_exciser(spec, exc_regs)

    assert_quantity_allclose(np.diff(excised_spec[50:60].flux),
                             np.diff(excised_spec[51:61].flux))
    assert_quantity_allclose(np.diff(excised_spec[80:90].flux),
                             np.diff(excised_spec[81:91].flux))


def test_extract_masked():
    wl = [1, 2, 3, 4]*u.nm
    flux = np.arange(4)*u.Jy
    mask = [False, False, True, True]

    masked_spec = Spectrum1D(spectral_axis=wl, flux=flux, mask=mask)
    region = SpectralRegion(1.5 * u.nm, 3.5 * u.nm)

    extracted = extract_region(masked_spec, region)

    assert np.all(extracted.mask == [False, True])
    assert np.all(extracted.flux.value == [1, 2])


def test_extract_multid_flux():
    flux = np.random.sample((10, 49)) * 100
    spec = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                      flux=flux * u.Jy)

    region = SpectralRegion(10 * u.nm, 20 * u.nm)
    extracted = extract_region(spec, region)

    assert extracted.shape == (10, 11)
    assert extracted[0,0].flux == spec[0,9].flux


def test_slab_multid_flux():
    flux = np.random.sample((10, 49)) * 100
    spec = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                      flux=flux * u.Jy)

    extracted = spectral_slab(spec, 10 * u.nm, 20 * u.nm)

    assert extracted.shape == (10, 11)
    assert extracted[0,0].flux == spec[0,9].flux
