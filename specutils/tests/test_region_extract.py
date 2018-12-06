import numpy as np
import pytest

import astropy.units as u
from astropy.nddata import StdDevUncertainty
from astropy.tests.helper import quantity_allclose

from ..spectra import Spectrum1D, SpectralRegion
from ..manipulation import extract_region
from ..manipulation.utils import linear_exciser
from .spectral_examples import simulated_spectra


def test_region_simple(simulated_spectra):
    np.random.seed(42)

    spectrum = simulated_spectra.s1_um_mJy_e1
    uncertainty = StdDevUncertainty(0.1*np.random.random(len(spectrum.flux))*u.mJy)
    spectrum.uncertainty = uncertainty

    region = SpectralRegion(0.6*u.um, 0.8*u.um)

    sub_spectrum = extract_region(spectrum, region)

    sub_spectrum_flux_expected = np.array(
            [1605.71612173, 1651.41650744, 2057.65798618, 2066.73502361, 1955.75832537,
             1670.52711471, 1491.10034446, 1637.08084112, 1471.28982259, 1299.19484483,
             1423.11195734, 1226.74494917, 1572.31888312, 1311.50503403, 1474.05051673,
             1335.39944397, 1420.61880528, 1433.18623759, 1290.26966668, 1605.67341284,
             1528.52281708, 1592.74392861, 1568.74162534, 1435.29407808, 1536.68040935,
             1157.33825995, 1136.12679394,  999.92394692, 1038.61546167, 1011.60297294])

    assert np.allclose(sub_spectrum.flux.value, sub_spectrum_flux_expected)


def test_region_ghz(simulated_spectra):
    spectrum = Spectrum1D(flux=simulated_spectra.s1_um_mJy_e1,
                          spectral_axis=simulated_spectra.s1_um_mJy_e1.frequency)

    region = SpectralRegion(374740.5725*u.GHz, 499654.09666667*u.GHz)

    sub_spectrum = extract_region(spectrum, region)

    sub_spectrum_flux_expected = np.array(
            [1605.71612173, 1651.41650744, 2057.65798618, 2066.73502361, 1955.75832537,
             1670.52711471, 1491.10034446, 1637.08084112, 1471.28982259, 1299.19484483,
             1423.11195734, 1226.74494917, 1572.31888312, 1311.50503403, 1474.05051673,
             1335.39944397, 1420.61880528, 1433.18623759, 1290.26966668, 1605.67341284,
             1528.52281708, 1592.74392861, 1568.74162534, 1435.29407808, 1536.68040935,
             1157.33825995, 1136.12679394,  999.92394692, 1038.61546167, 1011.60297294])

    assert np.allclose(sub_spectrum.flux.value, sub_spectrum_flux_expected)


def test_region_simple_check_ends(simulated_spectra):
    np.random.seed(42)

    spectrum = Spectrum1D(spectral_axis=np.linspace(1, 25, 25)*u.um, flux=np.random.random(25)*u.Jy)
    region = SpectralRegion(8*u.um, 15*u.um)

    sub_spectrum = extract_region(spectrum, region)

    assert sub_spectrum.spectral_axis.value[0] == 8
    assert sub_spectrum.spectral_axis.value[-1] == 15


def test_region_empty(simulated_spectra):
    np.random.seed(42)

    empty_spectrum = Spectrum1D(spectral_axis=[]*u.um, flux=[]*u.Jy)

    # Region past upper range of spectrum
    spectrum = Spectrum1D(spectral_axis=np.linspace(1, 25, 25)*u.um, flux=np.random.random(25)*u.Jy)
    region = SpectralRegion(28*u.um, 30*u.um)
    sub_spectrum = extract_region(spectrum, region)

    assert np.allclose(sub_spectrum.spectral_axis.value, empty_spectrum.spectral_axis.value)
    assert sub_spectrum.spectral_axis.unit == empty_spectrum.spectral_axis.unit

    assert np.allclose(sub_spectrum.flux.value, empty_spectrum.flux.value)
    assert sub_spectrum.flux.unit == empty_spectrum.flux.unit

    # Region below lower range of spectrum
    spectrum = Spectrum1D(spectral_axis=np.linspace(1, 25, 25)*u.um, flux=np.random.random(25)*u.Jy)
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
    with pytest.raises(Exception) as e_info:
        region = SpectralRegion(3*u.um, 3*u.um)


def test_region_two_sub(simulated_spectra):
    np.random.seed(42)

    spectrum = simulated_spectra.s1_um_mJy_e1
    uncertainty = StdDevUncertainty(0.1*np.random.random(len(spectrum.flux))*u.mJy)
    spectrum.uncertainty = uncertainty

    region = SpectralRegion([(0.6*u.um, 0.8*u.um), (0.86*u.um, 0.89*u.um)])

    sub_spectra = extract_region(spectrum, region)

    # Confirm the end points of the subspectra are correct
    assert np.allclose(sub_spectra[0].spectral_axis.value[[0, -1]],
                       np.array([0.6035353535353536, 0.793939393939394]))

    assert np.allclose(sub_spectra[1].spectral_axis.value[[0, -1]],
                       np.array([0.8661616161616162, 0.8858585858585859]))

    sub_spectrum_0_flux_expected = np.array(
            [1605.71612173, 1651.41650744, 2057.65798618, 2066.73502361, 1955.75832537,
             1670.52711471, 1491.10034446, 1637.08084112, 1471.28982259, 1299.19484483,
             1423.11195734, 1226.74494917, 1572.31888312, 1311.50503403, 1474.05051673,
             1335.39944397, 1420.61880528, 1433.18623759, 1290.26966668, 1605.67341284,
             1528.52281708, 1592.74392861, 1568.74162534, 1435.29407808, 1536.68040935,
             1157.33825995, 1136.12679394,  999.92394692, 1038.61546167, 1011.60297294])

    sub_spectrum_1_flux_expected = np.array(
            [1337.65312465, 1263.48914109, 1589.81797876, 1548.46068415])

    assert np.allclose(sub_spectra[0].flux.value, sub_spectrum_0_flux_expected)
    assert np.allclose(sub_spectra[1].flux.value, sub_spectrum_1_flux_expected)


def test_extract_region_pixels():
    spectrum = Spectrum1D(spectral_axis=np.linspace(4000, 10000, 25)*u.AA,
                          flux=np.arange(25)*u.Jy)
    region = SpectralRegion(10*u.pixel, 12*u.pixel)

    extracted = extract_region(spectrum, region)

    assert quantity_allclose(extracted.flux, [10, 11]*u.Jy)


def test_extract_region_mismatched_units():
    spectrum = Spectrum1D(spectral_axis=np.arange(25)*u.nm,
                          flux=np.arange(25)*u.Jy)

    region = SpectralRegion(100*u.AA, 119*u.AA)

    extracted = extract_region(spectrum, region)

    assert quantity_allclose(extracted.flux, [10, 11]*u.Jy)


def test_linear_excise_invert_from_spectrum():
    spec = Spectrum1D(flux=np.random.sample(100) * u.Jy,
                      spectral_axis=np.arange(100) * u.AA)
    inc_regs = SpectralRegion(80 * u.AA, 90 * u.AA) + \
               SpectralRegion(50 * u.AA, 60 * u.AA)
    exc_regs = inc_regs.invert_from_spectrum(spec)

    excised_spec = linear_exciser(spec, exc_regs)

    assert quantity_allclose(np.diff(excised_spec[50:60].flux),
                             np.diff(excised_spec[51:61].flux))
    assert quantity_allclose(np.diff(excised_spec[80:90].flux),
                             np.diff(excised_spec[81:91].flux))