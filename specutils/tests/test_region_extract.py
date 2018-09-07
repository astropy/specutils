import numpy as np

import astropy.units as u
from astropy.nddata import StdDevUncertainty

from ..spectra import Spectrum1D, SpectralRegion
from ..manipulation.extract_spectral_region import extract_region
from .spectral_examples import simulated_spectra


def test_region_simple(simulated_spectra):
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

    region = SpectralRegion(0.6*u.um, 0.8*u.um)

    sub_spectrum = extract_region(spectrum, region)

    sub_spectrum_flux_expected = np.array(
            [1605.71612173, 1651.41650744, 2057.65798618, 2066.73502361, 1955.75832537,
             1670.52711471, 1491.10034446, 1637.08084112, 1471.28982259, 1299.19484483,
             1423.11195734, 1226.74494917, 1572.31888312, 1311.50503403, 1474.05051673,
             1335.39944397, 1420.61880528, 1433.18623759, 1290.26966668, 1605.67341284,
             1528.52281708, 1592.74392861, 1568.74162534, 1435.29407808, 1536.68040935,
             1157.33825995, 1136.12679394,  999.92394692, 1038.61546167])

    assert np.allclose(sub_spectrum.flux.value, sub_spectrum_flux_expected)


def test_region_two_sub(simulated_spectra):
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

    region = SpectralRegion([(0.6*u.um, 0.8*u.um), (0.86*u.um, 0.89*u.um)])

    sub_spectrum = region.extract(spectrum)

    # TODO: Add in spectral_axis ck, but after PR #272 is in place

    sub_spectrum_flux_expected = np.array(
            [1605.71612173, 1651.41650744, 2057.65798618, 2066.73502361, 1955.75832537,
             1670.52711471, 1491.10034446, 1637.08084112, 1471.28982259, 1299.19484483,
             1423.11195734, 1226.74494917, 1572.31888312, 1311.50503403, 1474.05051673,
             1335.39944397, 1420.61880528, 1433.18623759, 1290.26966668, 1605.67341284,
             1528.52281708, 1592.74392861, 1568.74162534, 1435.29407808, 1536.68040935,
             1157.33825995, 1136.12679394,  999.92394692, 1038.61546167])

    assert np.allclose(sub_spectrum.flux.value, sub_spectrum_flux_expected)
