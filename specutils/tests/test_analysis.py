import pytest
import numpy as np

import astropy.units as u
from astropy.modeling import models
from astropy.nddata import StdDevUncertainty
from astropy.stats.funcs import gaussian_sigma_to_fwhm
from astropy.tests.helper import quantity_allclose

from ..spectra import Spectrum1D, SpectralRegion
from ..analysis import (line_flux, equivalent_width, snr, centroid,
                        gaussian_sigma_width, gaussian_fwhm, fwhm)
from ..tests.spectral_examples import simulated_spectra


def test_line_flux():

    np.random.seed(42)

    frequencies = np.linspace(1, 100, 10000) * u.GHz
    g = models.Gaussian1D(amplitude=1*u.Jy, mean=10*u.GHz, stddev=1*u.GHz)
    noise = np.random.normal(0., 0.01, frequencies.shape) * u.Jy
    flux = g(frequencies) + noise

    spectrum = Spectrum1D(spectral_axis=frequencies, flux=flux)

    result = line_flux(spectrum)

    assert result.unit.is_equivalent(u.erg / u.cm**2 / u.s)

    # Account for the fact that Astropy uses a different normalization of the
    # Gaussian where the integral is not 1
    expected = np.sqrt(2*np.pi) * u.GHz * u.Jy

    assert quantity_allclose(result, expected, atol=0.01*u.GHz*u.Jy)


def test_equivalent_width():

    np.random.seed(42)

    frequencies = np.linspace(1, 100, 10000) * u.GHz
    g = models.Gaussian1D(amplitude=1*u.Jy, mean=10*u.GHz, stddev=1*u.GHz)
    noise = np.random.normal(0., 0.01, frequencies.shape) * u.Jy
    flux = g(frequencies) + noise + 1*u.Jy

    spectrum = Spectrum1D(spectral_axis=frequencies, flux=flux)

    result = equivalent_width(spectrum)

    assert result.unit.is_equivalent(spectrum.spectral_axis_unit)

    # Since this is an emission line, we expect the equivalent width value to
    # be negative
    expected = -(np.sqrt(2*np.pi) * u.GHz)

    assert quantity_allclose(result, expected, atol=0.01*u.GHz)


def test_equivalent_width_regions():

    np.random.seed(42)

    frequencies = np.linspace(1, 100, 10000) * u.GHz
    g = models.Gaussian1D(amplitude=1*u.Jy, mean=10*u.GHz, stddev=1*u.GHz)
    noise = np.random.normal(0., 0.001, frequencies.shape) * u.Jy
    flux = g(frequencies) + noise + 1*u.Jy

    spec = Spectrum1D(spectral_axis=frequencies, flux=flux)
    cont_norm_spec = spec / np.median(spec.flux)
    result = equivalent_width(cont_norm_spec, regions=SpectralRegion(3*u.GHz, 97*u.GHz))

    expected = -(np.sqrt(2*np.pi) * u.GHz)

    assert quantity_allclose(result, expected, atol=0.01*u.GHz)


@pytest.mark.parametrize('continuum', [1*u.Jy, 2*u.Jy, 5*u.Jy])
def test_equivalent_width_continuum(continuum):

    np.random.seed(42)

    frequencies = np.linspace(1, 100, 10000) * u.GHz
    g = models.Gaussian1D(amplitude=1*u.Jy, mean=10*u.GHz, stddev=1*u.GHz)
    noise = np.random.normal(0., 0.01, frequencies.shape) * u.Jy
    flux = g(frequencies) + noise + continuum

    spectrum = Spectrum1D(spectral_axis=frequencies, flux=flux)

    result = equivalent_width(spectrum, continuum=continuum)

    assert result.unit.is_equivalent(spectrum.spectral_axis_unit)

    # Since this is an emission line, we expect the equivalent width value to
    # be negative
    expected = -(np.sqrt(2*np.pi) * u.GHz) / continuum.value

    assert quantity_allclose(result, expected, atol=0.01*u.GHz)


def test_equivalent_width_absorption():

    np.random.seed(42)

    frequencies = np.linspace(1, 100, 10000) * u.GHz
    amplitude = 0.5
    g = models.Gaussian1D(amplitude=amplitude*u.Jy, mean=10*u.GHz, stddev=1*u.GHz)
    continuum = 1*u.Jy
    noise = np.random.normal(0., 0.01, frequencies.shape) * u.Jy
    flux = continuum - g(frequencies) + noise

    spectrum = Spectrum1D(spectral_axis=frequencies, flux=flux)

    result = equivalent_width(spectrum)

    assert result.unit.is_equivalent(spectrum.spectral_axis_unit)

    # Since this is an absorption line, we expect the equivalent width value to
    # be positive
    expected = amplitude*np.sqrt(2*np.pi) * u.GHz

    assert quantity_allclose(result, expected, atol=0.01*u.GHz)


def test_snr(simulated_spectra):
    """
    Test the simple version of the spectral SNR.
    """

    np.random.seed(42)

    #
    #  Set up the data and add the uncertainty and calculate the expected SNR
    #

    spectrum = simulated_spectra.s1_um_mJy_e1
    uncertainty = StdDevUncertainty(0.1*np.random.random(len(spectrum.flux))*u.mJy)
    spectrum.uncertainty = uncertainty

    wavelengths = spectrum.spectral_axis
    flux = spectrum.flux

    spec_snr_expected = np.mean(flux / (uncertainty.array*uncertainty.unit))

    #
    # SNR of the whole spectrum
    #

    spec_snr = snr(spectrum)

    assert isinstance(spec_snr, u.Quantity)
    assert np.allclose(spec_snr.value, spec_snr_expected.value)


def test_snr_no_uncertainty(simulated_spectra):
    """
    Test the simple version of the spectral SNR.
    """

    #
    #  Set up the data and add the uncertainty and calculate the expected SNR
    #

    spectrum = simulated_spectra.s1_um_mJy_e1

    with pytest.raises(Exception) as e_info:
        _ = snr(spectrum)

def test_snr_multiple_flux(simulated_spectra):
    """
    Test the simple version of the spectral SNR, with multiple flux per single dispersion.
    """

    np.random.seed(42)

    #
    #  Set up the data and add the uncertainty and calculate the expected SNR
    #

    uncertainty = StdDevUncertainty(0.1*np.random.random((5, 10))*u.mJy)
    spec = Spectrum1D(spectral_axis=np.arange(10) * u.AA,
                      flux=np.random.sample((5, 10)) * u.Jy,
                      uncertainty=uncertainty)
    snr_spec = snr(spec)
    assert np.allclose(np.array(snr_spec), [18.20863867, 31.89475309, 14.51598119, 22.24603204, 32.01461421])

    uncertainty = StdDevUncertainty(0.1*np.random.random(10)*u.mJy)
    spec = Spectrum1D(spectral_axis=np.arange(10) * u.AA, flux=np.random.sample(10) * u.Jy, uncertainty=uncertainty)
    snr_spec = snr(spec)

    assert np.allclose(np.array(snr_spec), 31.325265361800415)


def test_snr_single_region(simulated_spectra):
    """
    Test the simple version of the spectral SNR over a region of the spectrum.
    """

    np.random.seed(42)

    region = SpectralRegion(0.52*u.um, 0.59*u.um)

    #
    #  Set up the data
    #

    spectrum = simulated_spectra.s1_um_mJy_e1
    uncertainty = StdDevUncertainty(0.1*np.random.random(len(spectrum.flux))*u.mJy)
    spectrum.uncertainty = uncertainty

    wavelengths = spectrum.spectral_axis
    flux = spectrum.flux

    # +1 because we want to include it in the calculation
    l = np.nonzero(wavelengths>region.lower)[0][0]
    r = np.nonzero(wavelengths<region.upper)[0][-1]+1

    spec_snr_expected = np.mean(flux[l:r] / (uncertainty.array[l:r]*uncertainty.unit))

    #
    # SNR of the whole spectrum
    #

    spec_snr = snr(spectrum, region)

    assert np.allclose(spec_snr.value, spec_snr_expected.value)


def test_snr_two_regions(simulated_spectra):
    """
    Test the simple version of the spectral SNR within two regions.
    """

    np.random.seed(42)

    #
    # Set the regions over which the SNR is calculated
    #

    regions = [SpectralRegion(0.52*u.um, 0.59*u.um), SpectralRegion(0.8*u.um, 0.9*u.um)]

    #
    #  Set up the data
    #

    spectrum = simulated_spectra.s1_um_mJy_e1
    uncertainty = StdDevUncertainty(0.1*np.random.random(len(spectrum.flux))*u.Jy)
    spectrum.uncertainty = uncertainty

    wavelengths = spectrum.spectral_axis
    flux = spectrum.flux

    spec_snr_expected = []
    for region in regions:

        l = np.nonzero(wavelengths >= region.lower)[0][0]
        r = np.nonzero(wavelengths <= region.upper)[0][-1]+1

        spec_snr_expected.append(np.mean(flux[l:r] / (uncertainty.array[l:r]*uncertainty.unit)))

    #
    # SNR of the whole spectrum
    #

    spec_snr = snr(spectrum, regions)

    assert np.allclose(spec_snr, spec_snr_expected)


def test_centroid(simulated_spectra):
    """
    Test the simple version of the spectral centroid.
    """

    np.random.seed(42)

    #
    #  Set up the data and add the uncertainty and calculate the expected SNR
    #

    spectrum = simulated_spectra.s1_um_mJy_e1
    uncertainty = StdDevUncertainty(0.1*np.random.random(len(spectrum.flux))*u.mJy)
    spectrum.uncertainty = uncertainty

    wavelengths = spectrum.spectral_axis
    flux = spectrum.flux

    spec_centroid_expected = np.sum(flux * wavelengths) / np.sum(flux)

    #
    # SNR of the whole spectrum
    #

    spec_centroid = centroid(spectrum, None)

    assert isinstance(spec_centroid, u.Quantity)
    assert np.allclose(spec_centroid.value, spec_centroid_expected.value)


def test_inverted_centroid(simulated_spectra):
    """
    Ensures the centroid calculation also works for *inverted* spectra - i.e.
    continuum-subtracted absorption lines.
    """
    spectrum = simulated_spectra.s1_um_mJy_e1
    spec_centroid_expected = (np.sum(spectrum.flux * spectrum.spectral_axis) /
                              np.sum(spectrum.flux))

    spectrum_inverted = Spectrum1D(spectral_axis=spectrum.spectral_axis,
                                   flux=-spectrum.flux)
    spec_centroid_inverted = centroid(spectrum_inverted, None)
    assert np.allclose(spec_centroid_inverted.value, spec_centroid_expected.value)


def test_centroid_multiple_flux(simulated_spectra):
    """
    Test the simple version of the spectral SNR, with multiple flux per single dispersion.
    """

    #
    #  Set up the data and add the uncertainty and calculate the expected SNR
    #

    np.random.seed(42)

    spec = Spectrum1D(spectral_axis=np.arange(10) * u.um,
                      flux=np.random.sample((5, 10)) * u.Jy)

    centroid_spec = centroid(spec, None)

    assert np.allclose(centroid_spec.value, np.array([4.46190995, 4.17223565, 4.37778249, 4.51595259, 4.7429066]))
    assert centroid_spec.unit == u.um


def test_gaussian_sigma_width():

    np.random.seed(42)

    # Create a (centered) gaussian spectrum for testing
    mean = 5
    frequencies = np.linspace(0, mean*2, 100) * u.GHz
    g1 = models.Gaussian1D(amplitude=5*u.Jy, mean=mean*u.GHz, stddev=0.8*u.GHz)

    spectrum = Spectrum1D(spectral_axis=frequencies, flux=g1(frequencies))

    result = gaussian_sigma_width(spectrum)

    assert quantity_allclose(result, g1.stddev, atol=0.01*u.GHz)


def test_gaussian_sigma_width_regions():

    np.random.seed(42)

    frequencies = np.linspace(0, 100, 10000) * u.GHz
    g1 = models.Gaussian1D(amplitude=5*u.Jy, mean=10*u.GHz, stddev=0.8*u.GHz)
    g2 = models.Gaussian1D(amplitude=5*u.Jy, mean=2*u.GHz, stddev=0.3*u.GHz)
    g3 = models.Gaussian1D(amplitude=5*u.Jy, mean=70*u.GHz, stddev=10*u.GHz)

    compound = g1 + g2 + g3
    spectrum = Spectrum1D(spectral_axis=frequencies, flux=compound(frequencies))

    region1 = SpectralRegion(5*u.GHz, 15*u.GHz)
    result1 = gaussian_sigma_width(spectrum, regions=region1)

    exp1 = g1.stddev
    assert quantity_allclose(result1, exp1, atol=0.25*exp1)

    region2 = SpectralRegion(1*u.GHz, 3*u.GHz)
    result2 = gaussian_sigma_width(spectrum, regions=region2)

    exp2 = g2.stddev
    assert quantity_allclose(result2, exp2, atol=0.25*exp2)

    region3 = SpectralRegion(40*u.GHz, 100*u.GHz)
    result3 = gaussian_sigma_width(spectrum, regions=region3)

    exp3 = g3.stddev
    assert quantity_allclose(result3, exp3, atol=0.25*exp3)

    # Test using a list of regions
    result_list = gaussian_sigma_width(spectrum, regions=[region1, region2, region3])
    for model, result in zip((g1, g2, g3), result_list):
        exp = model.stddev
        assert quantity_allclose(result, exp, atol=0.25*exp)


def test_gaussian_sigma_width_multi_spectrum():

    np.random.seed(42)

    frequencies = np.linspace(0, 100, 10000) * u.GHz
    g1 = models.Gaussian1D(amplitude=5*u.Jy, mean=50*u.GHz, stddev=0.8*u.GHz)
    g2 = models.Gaussian1D(amplitude=5*u.Jy, mean=50*u.GHz, stddev=5*u.GHz)
    g3 = models.Gaussian1D(amplitude=5*u.Jy, mean=50*u.GHz, stddev=10*u.GHz)

    flux = np.ndarray((3, len(frequencies))) * u.Jy

    flux[0] = g1(frequencies)
    flux[1] = g2(frequencies)
    flux[2] = g3(frequencies)

    spectra = Spectrum1D(spectral_axis=frequencies, flux=flux)

    results = gaussian_sigma_width(spectra)

    expected = (g1.stddev, g2.stddev, g3.stddev)
    for result, exp in zip(results, expected):
        assert quantity_allclose(result, exp, atol=0.25*exp)


def test_gaussian_fwhm():

    np.random.seed(42)

    # Create a (centered) gaussian spectrum for testing
    mean = 5
    frequencies = np.linspace(0, mean*2, 100) * u.GHz
    g1 = models.Gaussian1D(amplitude=5*u.Jy, mean=mean*u.GHz, stddev=0.8*u.GHz)

    spectrum = Spectrum1D(spectral_axis=frequencies, flux=g1(frequencies))

    result = gaussian_fwhm(spectrum)

    expected = g1.stddev * gaussian_sigma_to_fwhm
    assert quantity_allclose(result, expected, atol=0.01*u.GHz)


@pytest.mark.parametrize('mean', range(3,8))
def test_gaussian_fwhm_uncentered(mean):

    np.random.seed(42)

    # Create an uncentered gaussian spectrum for testing
    frequencies = np.linspace(0, 10, 1000) * u.GHz
    g1 = models.Gaussian1D(amplitude=5*u.Jy, mean=mean*u.GHz, stddev=0.8*u.GHz)

    spectrum = Spectrum1D(spectral_axis=frequencies, flux=g1(frequencies))

    result = gaussian_fwhm(spectrum)

    expected = g1.stddev * gaussian_sigma_to_fwhm
    assert quantity_allclose(result, expected, atol=0.05*u.GHz)


def test_fwhm():

    np.random.seed(42)

    # Create an (uncentered) spectrum for testing
    frequencies = np.linspace(0, 10, 1000) * u.GHz
    stddev = 0.8*u.GHz
    g1 = models.Gaussian1D(amplitude=5*u.Jy, mean=2*u.GHz, stddev=stddev)

    spectrum = Spectrum1D(spectral_axis=frequencies, flux=g1(frequencies))

    result = fwhm(spectrum)

    expected = stddev * gaussian_sigma_to_fwhm
    assert quantity_allclose(result, expected, atol=0.01*u.GHz)


def test_fwhm_multi_spectrum():

    np.random.seed(42)

    frequencies = np.linspace(0, 100, 10000) * u.GHz
    stddevs = [0.8, 5, 10]*u.GHz
    g1 = models.Gaussian1D(amplitude=5*u.Jy, mean=5*u.GHz, stddev=stddevs[0])
    g2 = models.Gaussian1D(amplitude=5*u.Jy, mean=50*u.GHz, stddev=stddevs[1])
    g3 = models.Gaussian1D(amplitude=5*u.Jy, mean=83*u.GHz, stddev=stddevs[2])

    flux = np.ndarray((3, len(frequencies))) * u.Jy

    flux[0] = g1(frequencies)
    flux[1] = g2(frequencies)
    flux[2] = g3(frequencies)

    spectra = Spectrum1D(spectral_axis=frequencies, flux=flux)

    results = fwhm(spectra)

    expected = stddevs * gaussian_sigma_to_fwhm
    assert quantity_allclose(results, expected, atol=0.01*u.GHz)
