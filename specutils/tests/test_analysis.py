import pytest
import numpy as np

import astropy.units as u
from astropy.modeling import models
from astropy.nddata import StdDevUncertainty
from astropy.stats.funcs import gaussian_sigma_to_fwhm
from astropy.tests.helper import quantity_allclose
from astropy.utils.exceptions import AstropyUserWarning

from ..spectra import Spectrum1D, SpectralRegion
from ..spectra.spectrum_collection import SpectrumCollection
from ..spectra.spectral_axis import SpectralAxis

from ..analysis import (line_flux, equivalent_width, snr, centroid,
                        gaussian_sigma_width, gaussian_fwhm, fwhm, moment,
                        snr_derived, fwzi, is_continuum_below_threshold)
from ..fitting import find_lines_threshold
from ..manipulation import snr_threshold, FluxConservingResampler


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


@pytest.mark.filterwarnings('ignore:invalid value encountered in true_divide:RuntimeWarning')
def test_line_flux_masked():

    np.random.seed(42)

    N = 100

    wavelengths = np.linspace(0.4, 1.05, N) * u.um

    g = models.Gaussian1D(amplitude=2000*u.mJy, mean=0.56*u.um, stddev=0.01*u.um)
    flux = g(wavelengths) + 1000 * u.mJy
    noise = 400 * np.random.random(flux.shape) * u.mJy
    flux += noise

    spectrum = Spectrum1D(spectral_axis=wavelengths, flux=flux)
    spectrum.uncertainty = StdDevUncertainty(noise)

    spectrum_masked = snr_threshold(spectrum, 10.)

    # Ensure we have at least 50% of the data being masked.
    assert len(np.where(spectrum_masked.mask)[0]) > N/2

    result = line_flux(spectrum_masked)

    assert result.unit.is_equivalent(u.Jy * u.um)

    assert quantity_allclose(result.value, 720.52992, atol=0.001)

    # With flux conserving resampler
    result = line_flux(spectrum_masked, mask_interpolation=FluxConservingResampler)
    assert quantity_allclose(result.value, 720.52992, atol=0.001)


def test_equivalent_width():

    np.random.seed(42)

    frequencies = np.linspace(100, 1, 10000) * u.GHz
    g = models.Gaussian1D(amplitude=1*u.Jy, mean=10*u.GHz, stddev=1*u.GHz)
    noise = np.random.normal(0., 0.01, frequencies.shape) * u.Jy
    flux = g(frequencies) + noise + 1*u.Jy

    spectrum = Spectrum1D(spectral_axis=frequencies, flux=flux)

    result = equivalent_width(spectrum)

    assert result.unit.is_equivalent(spectrum.wcs.unit, equivalencies=u.spectral())
    assert hasattr(result, 'uncertainty')

    # Since this is an emission line, we expect the equivalent width value to
    # be negative
    expected = -(np.sqrt(2*np.pi) * u.GHz)

    assert quantity_allclose(result, expected, atol=0.0025*u.GHz)


def test_equivalent_width_masked():

    np.random.seed(42)

    N = 100

    wavelengths = np.linspace(0.4, 1.05, N) * u.um

    g = models.Gaussian1D(amplitude=2000*u.mJy, mean=0.56*u.um, stddev=0.01*u.um)
    flux = g(wavelengths) + 1000 * u.mJy
    noise = 400 * np.random.random(flux.shape) * u.mJy
    flux += noise

    spectrum = Spectrum1D(spectral_axis=wavelengths, flux=flux)
    spectrum.uncertainty = StdDevUncertainty(noise)

    spectrum_masked = snr_threshold(spectrum, 10.)

    # Ensure we have at least 50% of the data being masked.
    assert len(np.where(spectrum_masked.mask)[0]) > N/2

    result = equivalent_width(spectrum_masked)

    assert result.unit.is_equivalent(spectrum.wcs.unit)

    assert quantity_allclose(result.value, -719.90618, atol=0.0005)

    # With flux conserving resampler
    result = equivalent_width(spectrum_masked,
                              mask_interpolation=FluxConservingResampler)
    assert quantity_allclose(result.value, -719.89962, atol=0.0005)


def test_equivalent_width_regions():

    np.random.seed(42)

    frequencies = np.linspace(100, 1, 10000) * u.GHz
    g = models.Gaussian1D(amplitude=1*u.Jy, mean=20*u.GHz, stddev=1*u.GHz)
    noise = np.random.normal(0., 0.001, frequencies.shape) * u.Jy
    flux = g(frequencies) + noise + 1*u.Jy

    spec = Spectrum1D(spectral_axis=frequencies, flux=flux)
    result = equivalent_width(spec, regions=SpectralRegion(60*u.GHz, 10*u.GHz))

    expected = -(np.sqrt(2*np.pi) * u.GHz)

    assert quantity_allclose(result, expected, atol=0.005*u.GHz)


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

    assert quantity_allclose(result, expected, atol=0.005*u.GHz)


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

    assert quantity_allclose(result, expected, atol=0.005*u.GHz)


@pytest.mark.parametrize('bin_specification', ["centers", "edges"])
def test_equivalent_width_bin_edges(bin_specification):
    """
    Test spectrum with bin specifications as centers or edges, and at modest sampling.
    """

    np.random.seed(42)

    wavelengths = np.linspace(500, 2500, 401) * u.AA
    spectral_axis = SpectralAxis(wavelengths, bin_specification=bin_specification)

    flunit = u.Unit('W m-2 AA-1')
    amplitude = 0.5
    stddev = 100*u.AA
    expected = amplitude*np.sqrt(2*np.pi) * stddev

    g = models.Gaussian1D(amplitude=amplitude*flunit, mean=1000*u.AA, stddev=stddev)
    continuum = 1 * flunit
    noise = np.random.normal(0., 0.01, spectral_axis.shape) * flunit
    flux = continuum - g(spectral_axis) + noise

    spectrum = Spectrum1D(spectral_axis=spectral_axis, flux=flux)

    result = equivalent_width(spectrum)

    assert result.unit.is_equivalent(spectrum.wcs.unit, equivalencies=u.spectral())
    assert quantity_allclose(result, expected, atol=0.5*u.AA)

    # With SpectralRegion; need higher tolerance as this is cutting more of the wings.
    result = equivalent_width(spectrum, regions=SpectralRegion(750*u.AA, 1500*u.AA))

    assert quantity_allclose(result, expected, atol=1*u.AA)


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

    flux = spectrum.flux

    spec_snr_expected = np.mean(flux / (uncertainty.array*uncertainty.unit))

    #
    # SNR of the whole spectrum
    #

    spec_snr = snr(spectrum)

    assert isinstance(spec_snr, u.Quantity)
    assert np.allclose(spec_snr.value, spec_snr_expected.value)


def test_snr_masked(simulated_spectra):
    """
    Test the simple version of the spectral SNR, with masked spectrum.
    """

    np.random.seed(42)

    spectrum = simulated_spectra.s1_um_mJy_e1_masked
    uncertainty = StdDevUncertainty(0.1*np.random.random(len(spectrum.flux))*u.mJy)
    spectrum.uncertainty = uncertainty

    uncertainty_array = uncertainty.array[~spectrum.mask]
    flux = spectrum.flux[~spectrum.mask]

    spec_snr_expected = np.mean(flux / (uncertainty_array * uncertainty.unit))

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

    with pytest.raises(Exception):
        snr(spectrum)


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
    l = np.nonzero(wavelengths > region.lower)[0][0]
    r = np.nonzero(wavelengths < region.upper)[0][-1]+1

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


def test_snr_masked_with_region(simulated_spectra):
    """
    Test the simple version of the spectral SNR over a region of the masked spectrum.
    """

    np.random.seed(42)

    spectrum = simulated_spectra.s1_um_mJy_e1_masked
    uncertainty = StdDevUncertainty(0.1 * np.random.random(len(spectrum.flux)) * u.mJy)
    spectrum.uncertainty = uncertainty

    wavelengths = spectrum.spectral_axis[~spectrum.mask]
    flux = spectrum.flux[~spectrum.mask]
    uncertainty_array = uncertainty.array[~spectrum.mask]

    region = SpectralRegion(0.52 * u.um, 0.59 * u.um)
    l = np.nonzero(wavelengths >= region.lower)[0][0]
    r = np.nonzero(wavelengths <= region.upper)[0][-1] + 1

    spec_snr_expected = np.mean(flux[l:r] / (uncertainty_array[l:r] * uncertainty.unit))

    spec_snr = snr(spectrum, region)

    assert np.allclose(spec_snr.value, spec_snr_expected.value)


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


def test_snr_derived_masked():
    np.random.seed(42)

    x = np.arange(1, 101) * u.um
    y = np.random.random(len(x))*u.Jy
    mask = (np.random.randn(x.shape[0])) > 0
    spectrum = Spectrum1D(spectral_axis=x, flux=y, mask=mask)

    assert np.allclose(snr_derived(spectrum), 2.08175408)

    sr = SpectralRegion(38*u.um, 48*u.um)
    assert np.allclose(snr_derived(spectrum, sr), 4.01610033)

    sr2 = SpectralRegion(48*u.um, 57*u.um)
    assert np.allclose(snr_derived(spectrum, [sr, sr2]), [4.01610033, 1.94906157])


@pytest.mark.parametrize("analytic", [True, False])
def test_centroid(simulated_spectra, analytic):
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

    spec_centroid = centroid(spectrum, analytic=analytic)

    assert isinstance(spec_centroid, u.Quantity)
    assert np.allclose(spec_centroid.value, spec_centroid_expected.value)
    assert hasattr(spec_centroid, 'uncertainty')
    # NOTE: value has not been scientifically validated
    if analytic:
        assert quantity_allclose(spec_centroid.uncertainty, 7.032035e-07*u.um, rtol=5e-5)
    else:
        assert quantity_allclose(spec_centroid.uncertainty, 6.916e-07*u.um, rtol=5e-3)


@pytest.mark.parametrize("analytic", [True, False])
def test_centroid_masked(simulated_spectra, analytic):
    """
    Test centroid with masked spectrum.
    """

    np.random.seed(42)

    # Same as in test for unmasked spectrum, but using
    # masked version of same spectrum.
    spectrum = simulated_spectra.s1_um_mJy_e1_masked
    uncertainty = StdDevUncertainty(0.1*np.random.random(len(spectrum.flux))*u.mJy)
    spectrum.uncertainty = uncertainty

    # Use masked flux and dispersion arrays to compute
    # the expected value for centroid.
    wavelengths = spectrum.spectral_axis[~spectrum.mask]
    flux = spectrum.flux[~spectrum.mask]

    spec_centroid_expected = np.sum(flux * wavelengths) / np.sum(flux)

    spec_centroid = centroid(spectrum, analytic=analytic)

    assert isinstance(spec_centroid, u.Quantity)
    assert np.allclose(spec_centroid.value, spec_centroid_expected.value)
    assert hasattr(spec_centroid, 'uncertainty')
    # NOTE: value has not been scientifically validated
    if analytic:
        assert quantity_allclose(spec_centroid.uncertainty, 1.87678219e-06*u.um, rtol=5e-5)
    else:
        assert quantity_allclose(spec_centroid.uncertainty, 1.92374917e-06*u.um, rtol=5e-5)


@pytest.mark.parametrize("analytic", [True, False])
def test_inverted_centroid(simulated_spectra, analytic):
    """
    Ensures the centroid calculation also works for *inverted* spectra - i.e.
    continuum-subtracted absorption lines.
    """
    spectrum = simulated_spectra.s1_um_mJy_e1
    uncertainty = StdDevUncertainty(0.1*np.random.random(len(spectrum.flux))*u.mJy)
    spectrum.uncertainty = uncertainty

    spec_centroid_expected = (np.sum(spectrum.flux * spectrum.spectral_axis) /
                              np.sum(spectrum.flux))

    spectrum_inverted = Spectrum1D(spectral_axis=spectrum.spectral_axis,
                                   flux=-spectrum.flux, uncertainty=uncertainty)
    spec_centroid_inverted = centroid(spectrum_inverted, analytic=analytic)
    assert np.allclose(spec_centroid_inverted.value, spec_centroid_expected.value)


@pytest.mark.parametrize("analytic", [True, False])
def test_inverted_centroid_masked(simulated_spectra, analytic):
    """
    Ensures the centroid calculation also works for *inverted* spectra with
    masked data - i.e. continuum-subtracted absorption lines.
    """
    spectrum = simulated_spectra.s1_um_mJy_e1_masked
    uncertainty = StdDevUncertainty(0.1*np.random.random(len(spectrum.flux))*u.mJy)
    spectrum.uncertainty = uncertainty

    spec_centroid_expected = (np.sum(spectrum.flux[~spectrum.mask] *
                                     spectrum.spectral_axis[~spectrum.mask]) /
                              np.sum(spectrum.flux[~spectrum.mask]))

    spectrum_inverted = Spectrum1D(spectral_axis=spectrum.spectral_axis,
                                   flux=-spectrum.flux,
                                   mask=spectrum.mask,
                                   uncertainty=uncertainty)

    spec_centroid_inverted = centroid(spectrum_inverted, analytic=analytic)
    assert np.allclose(spec_centroid_inverted.value, spec_centroid_expected.value)


@pytest.mark.parametrize("analytic", [True, False])
def test_centroid_multiple_flux(simulated_spectra, analytic):
    """
    Test the simple version of the spectral SNR, with multiple flux per single dispersion.
    """

    #
    #  Set up the data and add the uncertainty and calculate the expected SNR
    #

    np.random.seed(42)

    spec = Spectrum1D(spectral_axis=np.arange(1, 11) * u.um,
                      flux=np.random.sample((5, 10)) * u.Jy)
    uncertainty = StdDevUncertainty(0.1*np.random.random(spec.flux.shape)*u.mJy)
    spec.uncertainty = uncertainty

    centroid_spec = centroid(spec, analytic=analytic)
    print(centroid_spec.value)

    assert np.allclose(centroid_spec.value, np.array([5.46190995, 5.17223565, 5.37778249, 5.51595259, 5.7429066]))
    assert centroid_spec.unit == u.um
    assert hasattr(centroid_spec, 'uncertainty')
    assert len(centroid_spec.uncertainty) == 5
    # NOTE: value has not been scientifically validated
    if analytic:
        assert np.allclose(centroid_spec.uncertainty.value, np.array([1.14987628e-04, 1.49638658e-04, 1.02963584e-04, 1.21785134e-04, 9.47238087e-05]))
    else:
        assert np.allclose(centroid_spec.uncertainty.value, np.array([1.11040262e-04, 1.49617388e-04, 1.02730951e-04, 1.21124734e-04, 9.62905679e-05]))


@pytest.mark.parametrize("analytic", [True, False])
def test_gaussian_sigma_width(analytic):

    np.random.seed(42)

    # Create a (centered) gaussian spectrum for testing
    mean = 5
    frequencies = np.linspace(0, mean*2, 101)[1:] * u.GHz
    g1 = models.Gaussian1D(amplitude=5*u.Jy, mean=mean*u.GHz, stddev=0.8*u.GHz)

    spectrum = Spectrum1D(spectral_axis=frequencies, flux=g1(frequencies))
    uncertainty = StdDevUncertainty(0.1*np.random.random(len(spectrum.flux))*u.mJy)
    spectrum.uncertainty = uncertainty

    result = gaussian_sigma_width(spectrum, analytic=analytic)

    assert quantity_allclose(result, g1.stddev, atol=0.01*u.GHz)
    assert hasattr(result, 'uncertainty')
    # NOTE: value has not been scientifically validated!
    if analytic:
        assert quantity_allclose(result.uncertainty, 3.83737901e-05*u.GHz, rtol=5e-5)
    else:
        assert quantity_allclose(result.uncertainty, 3.59036951e-05*u.GHz, rtol=5e-5)


@pytest.mark.parametrize("analytic", [True, False])
def test_gaussian_sigma_width_masked(analytic):

    np.random.seed(42)

    # Create a (centered) gaussian masked spectrum for testing
    mean = 5
    frequencies = np.linspace(0, mean*2, 101)[1:] * u.GHz
    g1 = models.Gaussian1D(amplitude=5*u.Jy, mean=mean*u.GHz, stddev=0.8*u.GHz)
    uncertainty = StdDevUncertainty(0.1*np.random.random(len(frequencies))*u.Jy)

    mask = (np.random.randn(frequencies.shape[0]) + 1.) > 0

    spectrum = Spectrum1D(spectral_axis=frequencies, flux=g1(frequencies),
                          uncertainty=uncertainty, mask=mask)

    result = gaussian_sigma_width(spectrum, analytic=analytic)

    assert quantity_allclose(result, g1.stddev, atol=0.01*u.GHz)
    assert hasattr(result, 'uncertainty')
    # NOTE: value has not been scientifically validated!
    if analytic:
        assert quantity_allclose(result.uncertainty, 0.05245744*u.GHz, rtol=5e-5)
    else:
        assert quantity_allclose(result.uncertainty, 0.04954785*u.GHz, rtol=5e-5)


@pytest.mark.parametrize("analytic", [True, False])
def test_gaussian_sigma_width_regions(analytic):

    np.random.seed(42)

    frequencies = np.linspace(100, 0, 10000)[:-1] * u.GHz
    g1 = models.Gaussian1D(amplitude=5*u.Jy, mean=10*u.GHz, stddev=0.8*u.GHz)
    g2 = models.Gaussian1D(amplitude=5*u.Jy, mean=2*u.GHz, stddev=0.3*u.GHz)
    g3 = models.Gaussian1D(amplitude=5*u.Jy, mean=70*u.GHz, stddev=10*u.GHz)
    uncertainty = StdDevUncertainty(0.1*np.random.random(len(frequencies))*u.Jy)

    compound = g1 + g2 + g3
    spectrum = Spectrum1D(spectral_axis=frequencies, flux=compound(frequencies),
                          uncertainty=uncertainty)

    region1 = SpectralRegion(15*u.GHz, 5*u.GHz)
    result1 = gaussian_sigma_width(spectrum, regions=region1, analytic=analytic)

    exp1 = g1.stddev
    assert quantity_allclose(result1, exp1, atol=0.25*exp1)

    region2 = SpectralRegion(3*u.GHz, 1*u.GHz)
    result2 = gaussian_sigma_width(spectrum, regions=region2, analytic=analytic)

    exp2 = g2.stddev
    assert quantity_allclose(result2, exp2, atol=0.25*exp2)

    region3 = SpectralRegion(100*u.GHz, 40*u.GHz)
    result3 = gaussian_sigma_width(spectrum, regions=region3, analytic=analytic)

    exp3 = g3.stddev
    assert quantity_allclose(result3, exp3, atol=0.25*exp3)

    # Test using a list of regions
    result_list = gaussian_sigma_width(spectrum, regions=[region1, region2, region3],
                                       analytic=analytic)
    for model, result in zip((g1, g2, g3), result_list):
        exp = model.stddev
        assert quantity_allclose(result, exp, atol=0.25*exp)


@pytest.mark.parametrize("analytic", [True, False])
def test_gaussian_sigma_width_multi_spectrum(analytic):

    np.random.seed(42)

    frequencies = np.linspace(100, 0, 10000) * u.GHz
    g1 = models.Gaussian1D(amplitude=5*u.Jy, mean=50*u.GHz, stddev=0.8*u.GHz)
    g2 = models.Gaussian1D(amplitude=5*u.Jy, mean=50*u.GHz, stddev=5*u.GHz)
    g3 = models.Gaussian1D(amplitude=5*u.Jy, mean=50*u.GHz, stddev=10*u.GHz)

    flux = np.ndarray((3, len(frequencies))) * u.Jy

    flux[0] = g1(frequencies)
    flux[1] = g2(frequencies)
    flux[2] = g3(frequencies)

    # Add some noise so we don't have flux exactly = 0
    flux += 0.001*np.random.random(flux.shape)*u.Jy - 0.0005*u.Jy

    if analytic:
        spectra = Spectrum1D(spectral_axis=frequencies, flux=flux)
    else:
        uncertainty = StdDevUncertainty(0.1*np.random.random(flux.shape)*u.Jy)
        spectra = Spectrum1D(spectral_axis=frequencies, flux=flux, uncertainty=uncertainty)

    results = gaussian_sigma_width(spectra, analytic=analytic)

    expected = (g1.stddev, g2.stddev, g3.stddev)
    for result, exp in zip(results, expected):
        assert quantity_allclose(result, exp, atol=0.25*exp)


@pytest.mark.parametrize("analytic", [True, False])
def test_gaussian_fwhm(analytic):

    np.random.seed(42)

    # Create a (centered) gaussian spectrum for testing
    mean = 5
    frequencies = np.linspace(0, mean*2, 101)[1:] * u.GHz
    g1 = models.Gaussian1D(amplitude=5*u.Jy, mean=mean*u.GHz, stddev=0.8*u.GHz)

    spectrum = Spectrum1D(spectral_axis=frequencies, flux=g1(frequencies))
    uncertainty = StdDevUncertainty(0.1*np.random.random(len(spectrum.flux))*u.mJy)
    spectrum.uncertainty = uncertainty

    result = gaussian_fwhm(spectrum, analytic=analytic)

    expected = g1.stddev * gaussian_sigma_to_fwhm
    assert quantity_allclose(result, expected, atol=0.01*u.GHz)
    assert hasattr(result, 'uncertainty')
    # NOTE: value has not been scientifically validated!
    if analytic:
        assert quantity_allclose(result.uncertainty, 9.03633701e-05*u.GHz, rtol=5e-5)
    else:
        assert quantity_allclose(result.uncertainty, 8.45467409e-05*u.GHz, rtol=5e-5)


@pytest.mark.parametrize("analytic", [True, False])
def test_gaussian_fwhm_masked(analytic):

    np.random.seed(42)

    # Create a (centered) gaussian masked spectrum for testing
    mean = 5
    frequencies = np.linspace(0, mean*2, 100) * u.GHz
    g1 = models.Gaussian1D(amplitude=5*u.Jy, mean=mean*u.GHz, stddev=0.8*u.GHz)
    uncertainty = StdDevUncertainty(0.1*np.random.random(len(frequencies))*u.Jy)

    mask = (np.random.randn(frequencies.shape[0]) + 1.) > 0

    spectrum = Spectrum1D(spectral_axis=frequencies, flux=g1(frequencies),
                          uncertainty=uncertainty, mask=mask)

    result = gaussian_fwhm(spectrum, analytic=analytic)

    expected = g1.stddev * gaussian_sigma_to_fwhm
    assert quantity_allclose(result, expected, atol=0.01*u.GHz)
    assert hasattr(result, 'uncertainty')
    # NOTE: value has not been scientifically validated!
    if analytic:
        assert quantity_allclose(result.uncertainty, 0.12811291*u.GHz, rtol=5e-5)
    else:
        assert quantity_allclose(result.uncertainty, 0.12131201*u.GHz, rtol=5e-5)


@pytest.mark.parametrize('mean', range(3, 8))
def test_gaussian_fwhm_uncentered(mean):

    np.random.seed(42)

    # Create an uncentered gaussian spectrum for testing
    frequencies = np.linspace(0, 10, 1000) * u.GHz
    uncertainty=StdDevUncertainty(0.1*np.random.random(len(frequencies))*u.Jy)
    g1 = models.Gaussian1D(amplitude=5*u.Jy, mean=mean*u.GHz, stddev=0.8*u.GHz)

    spectrum = Spectrum1D(spectral_axis=frequencies, flux=g1(frequencies),
                          uncertainty=uncertainty)

    result = gaussian_fwhm(spectrum)

    expected = g1.stddev * gaussian_sigma_to_fwhm
    assert quantity_allclose(result, expected, atol=0.05*u.GHz)


def test_fwhm_masked():

    np.random.seed(42)

    # Create a masked (uncentered) spectrum for testing
    frequencies = np.linspace(0, 10, 1000) * u.GHz
    stddev = 0.8*u.GHz
    g1 = models.Gaussian1D(amplitude=5*u.Jy, mean=2*u.GHz, stddev=stddev)
    mask = (np.random.randn(frequencies.shape[0]) + 1.) > 0

    spectrum = Spectrum1D(spectral_axis=frequencies, flux=g1(frequencies),
                          mask=mask)

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
    flux = wavelengths.value * u.Jy  # highest point last.

    spectrum = Spectrum1D(spectral_axis=wavelengths, flux=flux)
    result = fwhm(spectrum)
    assert result == 5 * u.um

    # Flat spectrum
    wavelengths = np.linspace(1, 10, 100) * u.um
    flux = np.ones(wavelengths.shape) * u.Jy  # highest point last.

    spectrum = Spectrum1D(spectral_axis=wavelengths, flux=flux)
    result = fwhm(spectrum)
    assert result == 9 * u.um


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
    flux = wavelengths.value * u.Jy  # highest point last.

    spectrum = Spectrum1D(spectral_axis=wavelengths, flux=flux)
    result = fwhm(spectrum)
    assert result == 5*u.um

    # Flat spectrum
    wavelengths = np.linspace(1, 10, 100) * u.um
    flux = np.ones(wavelengths.shape) * u.Jy  # highest point last.

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


def test_fwzi_masked():
    np.random.seed(42)

    disp = np.linspace(0, 100, 100) * u.AA

    g = models.Gaussian1D(mean=np.mean(disp),
                          amplitude=1 * u.Jy,
                          stddev=10 * u.AA)
    flux = g(disp) + ((np.random.sample(disp.size) - 0.5) * 0.1) * u.Jy

    # Add mask. It is built such that about 50% of the data points
    # on and around the Gaussian peak are masked out (this was checked
    # with the debugger to examine in-memory data).
    uncertainty = StdDevUncertainty(0.1*np.random.random(len(disp))*u.Jy)
    mask = (np.random.randn(disp.shape[0]) - 0.5) > 0

    spec = Spectrum1D(spectral_axis=disp, flux=flux, uncertainty=uncertainty,
                      mask=mask)

    assert quantity_allclose(fwzi(spec), 35.9996284 * u.AA)


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

    assert is_continuum_below_threshold(spectrum, threshold=0.1*u.Jy) == True  # noqa

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

    assert is_continuum_below_threshold(spectrum, threshold=0.1*u.Jy) == True  # noqa

    # With mask, with uncertainty
    wavelengths = [300, 500, 1000] * u.nm
    data = [0.01, 1.029, 0.013] * u.Jy
    mask = np.array([False, True, False])
    uncertainty = StdDevUncertainty([1.01, 1.13, 1.1] * u.Jy)
    spectrum = Spectrum1D(spectral_axis=wavelengths, flux=data,
                          uncertainty=uncertainty, mask=mask)

    assert is_continuum_below_threshold(spectrum, threshold=0.1*u.Jy) == True  # noqa

    # With mask, with uncertainty -- should throw an exception
    wavelengths = [300, 500, 1000] * u.nm
    data = [10.03, 10.029, 10.033] * u.Jy
    mask = np.array([False, False, False])
    uncertainty = StdDevUncertainty([1.01, 1.13, 1.1] * u.Jy)
    spectrum = Spectrum1D(spectral_axis=wavelengths, flux=data,
                          uncertainty=uncertainty, mask=mask)

    with pytest.warns(AstropyUserWarning,
                      match=r'.*If you want to suppress this warning.*') as e_info:
        find_lines_threshold(spectrum, noise_factor=1)
    assert len(e_info) == 1


def test_moment():

    np.random.seed(42)

    frequencies = np.linspace(100, 1, 10000) * u.GHz
    g = models.Gaussian1D(amplitude=1*u.Jy, mean=10*u.GHz, stddev=1*u.GHz)
    noise = np.random.normal(0., 0.01, frequencies.shape) * u.Jy
    flux = g(frequencies) + noise

    spectrum = Spectrum1D(spectral_axis=frequencies, flux=flux)

    moment_0 = moment(spectrum, order=0)
    assert moment_0.unit.is_equivalent(u.Jy * u.GHz)
    assert quantity_allclose(moment_0, 2.5045*u.Jy*u.GHz, atol=0.001*u.Jy*u.GHz)

    moment_1 = moment(spectrum, order=1)
    assert moment_1.unit.is_equivalent(u.GHz)
    assert quantity_allclose(moment_1, 10.08*u.GHz, atol=0.01*u.GHz)

    moment_2 = moment(spectrum, order=2)
    assert moment_2.unit.is_equivalent(u.GHz**2)
    assert quantity_allclose(moment_2, 13.40*u.GHz**2, atol=0.01*u.GHz**2)

    moment_3 = moment(spectrum, order=3)
    assert moment_3.unit.is_equivalent(u.GHz**3)
    assert quantity_allclose(moment_3, 1233.78*u.GHz**3, atol=0.01*u.GHz**3)


def test_moment_cube():

    np.random.seed(42)

    frequencies = np.linspace(100, 1, 10000) * u.GHz
    g = models.Gaussian1D(amplitude=100*u.Jy, mean=50*u.GHz, stddev=1000*u.GHz)
    noise = np.random.normal(0., 0.01, frequencies.shape) * u.Jy
    flux = g(frequencies) + noise

    # do not use identical arrays in each spaxel. The purpose here is not to
    # check accuracy (already tested elsewhere), but dimensionality.

    flux_multid = np.broadcast_to(flux, [9, 10, flux.shape[0]]) * u.Jy

    spectrum = Spectrum1D(spectral_axis=frequencies, flux=flux_multid)

    moment_0 = moment(spectrum, order=0)

    assert moment_0.shape == (9, 10)
    assert moment_0.unit.is_equivalent(u.Jy*u.GHz)

    moment_1 = moment(spectrum, order=1)

    assert moment_1.shape == (9, 10)
    assert moment_1.unit.is_equivalent(u.GHz)
    assert quantity_allclose(moment_1, 50.50*u.GHz, atol=0.01*u.GHz)

    # select the last axis of the cube. Should be identical with
    # the default call above.
    moment_1_last = moment(spectrum, order=1, axis=2)

    assert moment_1_last.shape == moment_1.shape
    assert moment_1_last.unit.is_equivalent(moment_1.unit)
    assert quantity_allclose(moment_1, moment_1_last, rtol=1.E-5)

    # spatial 1st order - returns the dispersion
    moment_1 = moment(spectrum, order=1, axis=1)

    assert moment_1.shape == (9, 10000)
    assert moment_1.unit.is_equivalent(u.GHz)
    assert quantity_allclose(moment_1, frequencies, rtol=1.E-5)

    # higher order
    moment_2 = moment(spectrum, order=2)
    assert moment_2.shape == (9, 10)
    assert moment_2.unit.is_equivalent(u.GHz**2)
    assert quantity_allclose(moment_2, 816.648*u.GHz**2, atol=0.01*u.GHz**2)


def test_moment_cube_order_2():

    np.random.seed(42)

    frequencies = np.linspace(100, 1, 10000) * u.GHz
    g = models.Gaussian1D(amplitude=100*u.Jy, mean=50*u.GHz, stddev=1000*u.GHz)
    noise = np.random.normal(0., 0.01, frequencies.shape) * u.Jy
    flux = g(frequencies) + noise

    # NOTE: Does not work on non-square spatial-spatial slice.
    # use identical arrays in each spaxel. The purpose here is not to
    # check accuracy (already tested elsewhere), but dimensionality.

    flux_multid = np.broadcast_to(flux, [10, 10, flux.shape[0]]) * u.Jy

    spectrum = Spectrum1D(spectral_axis=frequencies, flux=flux_multid)

    # higher order
    moment_2 = moment(spectrum, order=2)

    assert moment_2.shape == (10, 10)
    assert moment_2.unit.is_equivalent(u.GHz**2)
    assert quantity_allclose(moment_2, 816.648*u.GHz**2, atol=0.01*u.GHz**2)

    # TODO: Remove test completely if it is agreed it is no longer needed
    # Github issue https://github.com/astropy/specutils/issues/930
    # mentions that removing ability to set axis makes this test invalid

    # spatial higher order (what's the meaning of this?)
    moment_2 = moment(spectrum, order=2, axis=1)

    assert moment_2.shape == (10, 10000)
    assert moment_2.unit.is_equivalent(u.GHz**2)
    # check assorted values.
    assert quantity_allclose(moment_2[0][0], 8.078e-28*u.GHz**2, rtol=0.01)
    assert quantity_allclose(moment_2[1][0], 8.078e-28*u.GHz**2, rtol=0.01)
    assert quantity_allclose(moment_2[0][3], 2.019e-28*u.GHz**2, rtol=0.01)


def test_moment_cube_order_1_to_6():

    np.random.seed(42)

    frequencies = np.linspace(100, 1, 10000) * u.Hz
    sig = 5
    mean = 50
    amp = 1./(sig*np.sqrt(2*np.pi))
    g = models.Gaussian1D(amplitude=amp, mean=mean*u.Hz, stddev=sig*u.Hz)
    flux = g(frequencies)

    flux_multid = np.broadcast_to(flux, [9, 10, flux.shape[0]]) * u.Jy

    spectrum = Spectrum1D(spectral_axis=frequencies, flux=flux_multid)

    moment_1 = moment(spectrum, order=1)

    assert moment_1.shape == (9, 10)
    assert moment_1.unit.is_equivalent(u.Hz)
    assert quantity_allclose(moment_1, 50.0*u.Hz, atol=0.01*u.Hz)

    moment_2 = moment(spectrum, order=2)
    assert moment_2.shape == (9, 10)
    assert moment_2.unit.is_equivalent(u.Hz**2)
    assert quantity_allclose(moment_2, 25.0*u.Hz**2, atol=0.01*u.Hz**2)

    moment_3 = moment(spectrum, order=3)
    assert moment_3.shape == (9, 10)
    assert moment_3.unit.is_equivalent(u.Hz**3)
    assert quantity_allclose(moment_3, 0.0*u.Hz**3, atol=0.01*u.Hz**3)

    moment_4 = moment(spectrum, order=4)
    assert moment_4.shape == (9, 10)
    assert moment_4.unit.is_equivalent(u.Hz**4)
    assert quantity_allclose(moment_4, 1875.0*u.Hz**4, atol=0.01*u.Hz**4)

    moment_5 = moment(spectrum, order=5)
    assert moment_5.shape == (9, 10)
    assert moment_5.unit.is_equivalent(u.Hz**5)
    assert quantity_allclose(moment_5, 0.0*u.Hz**5, atol=0.01*u.Hz**5)

    moment_6 = moment(spectrum, order=6)
    assert moment_6.shape == (9, 10)
    assert moment_6.unit.is_equivalent(u.Hz**6)
    assert quantity_allclose(moment_6, 234375.0*u.Hz**6, atol=0.01*u.Hz**6)


@pytest.mark.filterwarnings('ignore:Not all spectra have associated uncertainties of the same type')
@pytest.mark.filterwarnings('ignore:Not all spectra have associated masks')
def test_moment_collection():

    np.random.seed(42)

    frequencies = np.linspace(100, 1, 10000) * u.GHz
    g = models.Gaussian1D(amplitude=1*u.Jy, mean=10*u.GHz, stddev=1*u.GHz)
    noise = np.random.normal(0., 0.01, frequencies.shape) * u.Jy
    flux = g(frequencies) + noise
    s1 = Spectrum1D(spectral_axis=frequencies, flux=flux)

    frequencies = np.linspace(200, 2, 10000) * u.GHz
    g = models.Gaussian1D(amplitude=2*u.Jy, mean=20*u.GHz, stddev=2*u.GHz)
    noise = np.random.normal(0., 0.02, frequencies.shape) * u.Jy
    flux = g(frequencies) + noise
    s2 = Spectrum1D(spectral_axis=frequencies, flux=flux)

    collection = SpectrumCollection.from_spectra([s1, s2])

    # Compare moments derived from collection, with moments
    # derived from individual members.
    moment_0 = moment(collection, order=0)
    moment_0_s1 = moment(s1, order=0)
    moment_0_s2 = moment(s2, order=0)
    assert quantity_allclose(moment_0[0], moment_0_s1, rtol=1.E-5)
    assert quantity_allclose(moment_0[1], moment_0_s2, rtol=1.E-5)

    moment_1 = moment(collection, order=1)
    moment_1_s1 = moment(s1, order=1)
    moment_1_s2 = moment(s2, order=1)
    assert quantity_allclose(moment_1[0], moment_1_s1, rtol=1.E-5)
    assert quantity_allclose(moment_1[1], moment_1_s2, rtol=1.E-5)

    moment_2 = moment(collection, order=2)
    moment_2_s1 = moment(s1, order=2)
    moment_2_s2 = moment(s2, order=2)
    assert quantity_allclose(moment_2[0], moment_2_s1, rtol=1.E-5)
    assert quantity_allclose(moment_2[1], moment_2_s2, rtol=1.E-5)
