import astropy.units as u
import astropy.wcs as fitswcs
import gwcs
import numpy as np
import pytest

from ..tests.spectral_examples import simulated_spectra
from ..spectra.spectrum1d import Spectrum1D
from astropy.nddata import StdDevUncertainty
import astropy.units as u
from ..analysis import equivalent_width, snr
from ..utils import SpectralRegion


def test_equivalent_width():
    spec = Spectrum1D(spectral_axis=np.arange(50),
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
        0.19552929, -0.7862626 ,  0.51592792, -0.95650517, -1.26917689]))

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

    l = np.nonzero(wavelengths>region.lower)[0][0]
    r = np.nonzero(wavelengths<region.upper)[0][-1]

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

        l = np.nonzero(wavelengths>region.lower)[0][0]
        r = np.nonzero(wavelengths<region.upper)[0][-1]

        spec_snr_expected.append(np.mean(flux[l:r] / (uncertainty.array[l:r]*uncertainty.unit)))

    #
    # SNR of the whole spectrum
    #

    spec_snr = snr(spectrum, regions)

    assert np.allclose(spec_snr, spec_snr_expected)


def test_snr_single_region_with_noise_region(simulated_spectra):
    """
    Test the simple version of the spectral SNR over a region of the spectrum.
    """

    np.random.seed(42)

    region = SpectralRegion(0.52*u.um, 0.59*u.um)
    noise_region = SpectralRegion(0.40*u.um, 0.45*u.um)
    
    #
    #  Set up the data
    #

    spectrum = simulated_spectra.s1_um_mJy_e1
    uncertainty = StdDevUncertainty(0.1*np.random.random(len(spectrum.flux))*u.mJy)
    spectrum.uncertainty = uncertainty

    wavelengths = spectrum.spectral_axis
    flux = spectrum.flux

    l = np.nonzero(wavelengths>=region.lower)[0][0]
    r = np.nonzero(wavelengths<region.upper)[0][-1]

    noise_l = np.nonzero(wavelengths>=noise_region.lower)[0][0]
    noise_r = np.nonzero(wavelengths<noise_region.upper)[0][-1]

    spec_snr_expected = np.mean(flux[l:r] / np.std(flux[noise_l:noise_r]))

    #
    # SNR of the whole spectrum
    #

    spec_snr = snr(spectrum, region, noise_region=noise_region)

    assert np.allclose(spec_snr.value, spec_snr_expected.value)


