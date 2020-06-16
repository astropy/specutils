import pytest
import numpy as np

import astropy.units as u
from astropy.modeling import models
from astropy.nddata import StdDevUncertainty
from astropy.stats.funcs import gaussian_sigma_to_fwhm
from astropy.tests.helper import quantity_allclose
from astropy.utils.exceptions import AstropyUserWarning

from ..spectra import Spectrum1D, SpectralRegion
from ..analysis import (line_flux, equivalent_width, snr, centroid,
                        gaussian_sigma_width, gaussian_fwhm, fwhm,
                        snr_derived, fwzi, is_continuum_below_threshold)
from ..fitting import find_lines_threshold
from ..tests.spectral_examples import simulated_spectra


def test_line_flux():

    np.random.seed(42)

    frequencies = np.linspace(100, 1, 10000) * u.GHz
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


def test_line_flux_uncertainty():
    np.random.seed(42)

    spec = Spectrum1D(spectral_axis=np.arange(10) * u.AA,
                      flux=np.random.sample(10) * u.Jy,
                      uncertainty=StdDevUncertainty(np.random.sample(10) * 0.01))

    lf = line_flux(spec)

    assert quantity_allclose(lf, u.Quantity(5.20136736, u.Jy * u.AA))
    assert quantity_allclose(lf.uncertainty, u.Quantity(0.01544415, u.Jy * u.AA))


def test_equivalent_width():

    np.random.seed(42)

    frequencies = np.linspace(100, 1, 10000) * u.GHz
    g = models.Gaussian1D(amplitude=1*u.Jy, mean=10*u.GHz, stddev=1*u.GHz)
    noise = np.random.normal(0., 0.01, frequencies.shape) * u.Jy
    flux = g(frequencies) + noise + 1*u.Jy

    spectrum = Spectrum1D(spectral_axis=frequencies, flux=flux)

    result = equivalent_width(spectrum)

    assert result.unit.is_equivalent(spectrum.wcs.unit, equivalencies=u.spectral())

    # Since this is an emission line, we expect the equivalent width value to
    # be negative
    expected = -(np.sqrt(2*np.pi) * u.GHz)

    assert quantity_allclose(result, expected, atol=0.01*u.GHz)


def test_equivalent_width_regions():

    np.random.seed(42)

    frequencies = np.linspace(100, 1, 10000) * u.GHz
    g = models.Gaussian1D(amplitude=1*u.Jy, mean=10*u.GHz, stddev=1*u.GHz)
    noise = np.random.normal(0., 0.001, frequencies.shape) * u.Jy
    flux = g(frequencies) + noise + 1*u.Jy

    spec = Spectrum1D(spectral_axis=frequencies, flux=flux)
    cont_norm_spec = spec / np.median(spec.flux)
    result = equivalent_width(cont_norm_spec, regions=SpectralRegion(97*u.GHz, 3*u.GHz))

    expected = -(np.sqrt(2*np.pi) * u.GHz)

    assert quantity_allclose(result, expected, atol=0.02*u.GHz)


@pytest.mark.parametrize('continuum', [1*u.Jy, 2*u.Jy, 5*u.Jy])
def test_equivalent_width_continuum(continuum):

    np.random.seed(42)

    frequencies = np.linspace(100, 1, 10000) * u.GHz
    g = models.Gaussian1D(amplitude=1*u.Jy, mean=10*u.GHz, stddev=1*u.GHz)
    noise = np.random.normal(0., 0.01, frequencies.shape) * u.Jy
    flux = g(frequencies) + noise + continuum

    spectrum = Spectrum1D(spectral_axis=frequencies, flux=flux)

    result = equivalent_width(spectrum, continuum=continuum)

    assert result.unit.is_equivalent(spectrum.wcs.unit, equivalencies=u.spectral())

    # Since this is an emission line, we expect the equivalent width value to
    # be negative
    expected = -(np.sqrt(2*np.pi) * u.GHz) / continuum.value

    assert quantity_allclose(result, expected, atol=0.01*u.GHz)


def test_equivalent_width_absorption():

    np.random.seed(42)

    frequencies = np.linspace(100, 1, 10000) * u.GHz
    amplitude = 0.5
    g = models.Gaussian1D(amplitude=amplitude*u.Jy, mean=10*u.GHz, stddev=1*u.GHz)
    continuum = 1*u.Jy
    noise = np.random.normal(0., 0.01, frequencies.shape) * u.Jy
    flux = continuum - g(frequencies) + noise

    spectrum = Spectrum1D(spectral_axis=frequencies, flux=flux)

    result = equivalent_width(spectrum)

    assert result.unit.is_equivalent(spectrum.wcs.unit, equivalencies=u.spectral())

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


def test_snr_derived():
    np.random.seed(42)

    x = np.arange(1, 101) * u.um
    y = np.random.random(len(x))*u.Jy
    spectrum = Spectrum1D(spectral_axis=x, flux=y)

    assert np.allclose(snr_derived(spectrum), 1.604666860424951)

    sr = SpectralRegion(38*u.um, 48*u.um)
    assert np.allclose(snr_derived(spectrum, sr), 2.330463630828406)

    sr2 = SpectralRegion(48*u.um, 57*u.um)
    assert np.allclose(snr_derived(spectrum, [sr, sr2]), [2.330463630828406, 2.9673559890209305])


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

    frequencies = np.linspace(100, 0, 10000) * u.GHz
    g1 = models.Gaussian1D(amplitude=5*u.Jy, mean=10*u.GHz, stddev=0.8*u.GHz)
    g2 = models.Gaussian1D(amplitude=5*u.Jy, mean=2*u.GHz, stddev=0.3*u.GHz)
    g3 = models.Gaussian1D(amplitude=5*u.Jy, mean=70*u.GHz, stddev=10*u.GHz)

    compound = g1 + g2 + g3
    spectrum = Spectrum1D(spectral_axis=frequencies, flux=compound(frequencies))

    region1 = SpectralRegion(15*u.GHz, 5*u.GHz)
    result1 = gaussian_sigma_width(spectrum, regions=region1)

    exp1 = g1.stddev
    assert quantity_allclose(result1, exp1, atol=0.25*exp1)

    region2 = SpectralRegion(3*u.GHz, 1*u.GHz)
    result2 = gaussian_sigma_width(spectrum, regions=region2)

    exp2 = g2.stddev
    assert quantity_allclose(result2, exp2, atol=0.25*exp2)

    region3 = SpectralRegion(100*u.GHz, 40*u.GHz)
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

    frequencies = np.linspace(100, 0, 10000) * u.GHz
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

    # Highest point at the first point
    wavelengths = np.linspace(1, 10, 100) * u.um
    flux = (1.0 / wavelengths.value) * u.Jy  # highest point first.

    spectrum = Spectrum1D(spectral_axis=wavelengths, flux=flux)
    result = fwhm(spectrum)
    # Note that this makes a little more sense than the previous version;
    # since the maximum value occurs at wavelength=1, and the half-value of
    # flux (0.5) occurs at exactly wavelength=2, the result should be
    # exactly 1 (2 - 1).
    assert result == 1.0 * u.um

    # Test the interpolation used in FWHM for wavelength values that are not
    # on the grid
    wavelengths = np.linspace(1, 10, 31) * u.um
    flux = (1.0 / wavelengths.value) * u.Jy  # highest point first.

    spectrum = Spectrum1D(spectral_axis=wavelengths, flux=flux)
    result = fwhm(spectrum)

    assert quantity_allclose(result, 1.01 * u.um)

    # Highest point at the last point
    wavelengths = np.linspace(1, 10, 100) * u.um
    flux = wavelengths.value*u.Jy # highest point last.

    spectrum = Spectrum1D(spectral_axis=wavelengths, flux=flux)
    result = fwhm(spectrum)
    assert result == 5*u.um

    # Flat spectrum
    wavelengths = np.linspace(1, 10, 100) * u.um
    flux = np.ones(wavelengths.shape)*u.Jy # highest point last.

    spectrum = Spectrum1D(spectral_axis=wavelengths, flux=flux)
    result = fwhm(spectrum)
    assert result == 9*u.um


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


def test_fwzi():
    np.random.seed(42)

    disp = np.linspace(0, 100, 1000) * u.AA

    g = models.Gaussian1D(mean=np.mean(disp),
                          amplitude=1 * u.Jy,
                          stddev=2 * u.AA)

    flux = g(disp)
    flux_with_noise = g(disp) + ((np.random.sample(disp.size) - 0.5) * 0.1) * u.Jy

    spec = Spectrum1D(spectral_axis=disp, flux=flux)
    spec_with_noise = Spectrum1D(spectral_axis=disp, flux=flux_with_noise)

    assert quantity_allclose(fwzi(spec), 226.89732509 * u.AA)
    assert quantity_allclose(fwzi(spec_with_noise), 106.99998944 * u.AA)


def test_fwzi_multi_spectrum():
    np.random.seed(42)

    disp = np.linspace(0, 100, 1000) * u.AA

    amplitudes = [0.1, 1, 10] * u.Jy
    means = [25, 50, 75] * u.AA
    stddevs = [1, 5, 10] * u.AA
    params = list(zip(amplitudes, means, stddevs))

    flux = np.zeros(shape=(3, len(disp)))

    for i in range(3):
        flux[i] = models.Gaussian1D(*params[i])(disp)

    spec = Spectrum1D(spectral_axis=disp, flux=flux * u.Jy)

    expected = [113.51706001 * u.AA, 567.21252727 * u.AA, 499.5024546 * u.AA]

    assert quantity_allclose(fwzi(spec), expected)


def test_is_continuum_below_threshold():

    # No mask, no uncertainty
    wavelengths = [300, 500, 1000] * u.nm
    data = [0.001, -0.003, 0.003] * u.Jy
    spectrum = Spectrum1D(spectral_axis=wavelengths, flux=data)
    assert True == is_continuum_below_threshold(spectrum, threshold=0.1*u.Jy)

#    # No mask, no uncertainty, threshold is float
#    wavelengths = [300, 500, 1000] * u.nm
#    data = [0.0081, 0.0043, 0.0072] * u.Jy
#    spectrum = Spectrum1D(spectral_axis=wavelengths, flux=data)
#    assert True == is_continuum_below_threshold(spectrum, threshold=0.1)

    # No mask, with uncertainty
    wavelengths = [300, 500, 1000] * u.nm
    data = [0.03, 0.029, 0.031] * u.Jy
    uncertainty = StdDevUncertainty([1.01, 1.03, 1.01] * u.Jy)
    spectrum = Spectrum1D(spectral_axis=wavelengths, flux=data,
                          uncertainty=uncertainty)

    assert True == is_continuum_below_threshold(spectrum, threshold=0.1*u.Jy)

    # With mask, with uncertainty
    wavelengths = [300, 500, 1000] * u.nm
    data = [0.01, 1.029, 0.013] * u.Jy
    mask = np.array([False, True, False])
    uncertainty = StdDevUncertainty([1.01, 1.13, 1.1] * u.Jy)
    spectrum = Spectrum1D(spectral_axis=wavelengths, flux=data,
                          uncertainty=uncertainty, mask=mask)

    assert True == is_continuum_below_threshold(spectrum, threshold=0.1*u.Jy)

    # With mask, with uncertainty -- should throw an exception
    wavelengths = [300, 500, 1000] * u.nm
    data = [10.03, 10.029, 10.033] * u.Jy
    mask = np.array([False, False, False])
    uncertainty = StdDevUncertainty([1.01, 1.13, 1.1] * u.Jy)
    spectrum = Spectrum1D(spectral_axis=wavelengths, flux=data,
                          uncertainty=uncertainty, mask=mask)

    print('spectrum has flux {}'.format(spectrum.flux))
    with pytest.warns(AstropyUserWarning) as e_info:
        find_lines_threshold(spectrum, noise_factor=1)
        assert len(e_info)==1 and 'if you want to suppress this warning' in e_info[0].message.args[0].lower()
