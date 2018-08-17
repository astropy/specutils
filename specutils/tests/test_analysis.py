import pytest
import numpy as np

import astropy.units as u
from astropy.units import quantity
from astropy.modeling import models
from astropy.nddata import StdDevUncertainty

from ..spectra import Spectrum1D, SpectralRegion
from ..analysis import equivalent_width, snr, centroid, sigma_full_width
from ..manipulation import noise_region_uncertainty
from ..tests.spectral_examples import simulated_spectra


def test_equivalent_width():
    spec = Spectrum1D(spectral_axis=np.arange(50) * u.AA,
                      flux=np.array([
        0.09731049,  0.04101204,  0.09726903,  0.72524453, -0.28105335,
       -0.93772961, -1.9828759 ,  0.38752423,  0.86006845, -0.08198352,
       -0.08303639,  0.18421212, -0.50724803, -0.09625829, -1.9318252 ,
        1.07973092, -0.3189966 , -1.52045995,  1.95926732,  1.71674612,
        0.28450979,  0.37737352, -1.16525665,  0.29277855, -0.37458935,
       -1.31719473, -0.31894975, -0.51095169, -0.45959643,  0.77837891,
        0.91153499,  0.13612405,  0.63433898, -0.91986964, -0.4546604 ,
       -1.09797558, -1.83641516, -0.94179757, -1.33989398, -0.66452815,
       -0.71835507, -1.39939311,  0.50070437, -1.03926682,  0.58481419,
        0.19552929, -0.7862626 ,  0.51592792, -0.95650517, -1.26917689]) * u.Jy)

    ew = equivalent_width(spec)

    assert isinstance(ew, u.Quantity)
    assert np.allclose(ew.value, 6.8278704893358)

    ew = equivalent_width(spec[10:20])

    assert np.allclose(ew.value, 15.37809622)


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


def test_sigma_full_width():

    np.random.seed(42)

    # Create a (centered) gaussian spectrum for testing
    mean = 5
    frequencies = np.linspace(1, mean*2, 100) * u.GHz
    g1 = models.Gaussian1D(amplitude=5*u.Jy, mean=mean*u.GHz, stddev=0.8*u.GHz)

    spectrum = Spectrum1D(spectral_axis=frequencies, flux=g1(frequencies))

    result = sigma_full_width(spectrum)

    assert quantity.isclose(result, g1.stddev*2, atol=0.01*u.GHz)
