import numpy as np
from astropy import convolution
import astropy.units as u
from ..spectra.spectrum1d import Spectrum1D

from ..processing.smoothing import (box_smooth, gaussian_smooth,
                                    mexicanhat_smooth, trapezoid_smooth)


def test_smooth_box():

    #  Create the spectrum
    spec1 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                       flux=np.random.sample(49) * 100)

    flux = spec1.flux
    width = 3

    # Create the flux_smoothed which is what we want to compare to
    box_kernel = convolution.Gaussian1D(width)
    flux_smoothed = convolution.convolve(flux, box_kernel)

    # Test not-in-place smoothing, defualt is not in-place
    spec1_smoothed = box_smooth(spec1, width)
    assert np.allclose(spec1_smoothed, flux_smoothed)

    # Test not-in-place smoothing
    spec1_smoothed = box_smooth(spec1, width, inplace=False)
    assert np.allclose(spec1_smoothed, flux_smoothed)

    # Test in-place smoothing
    spec1_smoothed = box_smooth(spec1, width, inplace=True)
    assert np.allclose(spec1_smoothed, flux_smoothed)


def test_smooth_gaussian():

    #  Create the spectrum
    spec1 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                       flux=np.random.sample(49) * 100)

    flux = spec1.flux
    stddev = 3

    # Create the flux_smoothed which is what we want to compare to
    gaussian_kernel = convolution.Gaussian1D(stddev)
    flux_smoothed = convolution.convolve(flux, gaussian_kernel)

    # Test not-in-place smoothing, defualt is not in-place
    spec1_smoothed = gaussian_smooth(spec1, stddev)
    assert np.allclose(spec1_smoothed, flux_smoothed)

    # Test not-in-place smoothing
    spec1_smoothed = gaussian_smooth(spec1, stddev, inplace=False)
    assert np.allclose(spec1_smoothed, flux_smoothed)

    # Test in-place smoothing
    spec1_smoothed = gaussian_smooth(spec1, stddev, inplace=True)
    assert np.allclose(spec1_smoothed, flux_smoothed)


def test_smooth_mexicanhat():

    #  Create the spectrum
    spec1 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                       flux=np.random.sample(49) * 100)

    flux = spec1.flux
    stddev = 3

    # Create the flux_smoothed which is what we want to compare to
    mexicanhat_kernel = convolution.MexicanHat1D(stddev)
    flux_smoothed = convolution.convolve(flux, mexicanhat_kernel)

    # Test not-in-place smoothing, defualt is not in-place
    spec1_smoothed = mexicanhat_smooth(spec1, stddev)
    assert np.allclose(spec1_smoothed, flux_smoothed)

    # Test not-in-place smoothing
    spec1_smoothed = mexicanhat_smooth(spec1, stddev, inplace=False)
    assert np.allclose(spec1_smoothed, flux_smoothed)

    # Test in-place smoothing
    spec1_smoothed = mexicanhat_smooth(spec1, stddev, inplace=True)
    assert np.allclose(spec1_smoothed, flux_smoothed)


def test_smooth_trapezoid():

    #  Create the spectrum
    spec1 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                       flux=np.random.sample(49) * 100)

    flux = spec1.flux
    stddev = 3

    # Create the flux_smoothed which is what we want to compare to
    trapezoid_kernel = convolution.Trapezoid1D(stddev)
    flux_smoothed = convolution.convolve(flux, trapezoid_kernel)

    # Test not-in-place smoothing, defualt is not in-place
    spec1_smoothed = trapezoid_smooth(spec1, stddev)
    assert np.allclose(spec1_smoothed, flux_smoothed)

    # Test not-in-place smoothing
    spec1_smoothed = trapezoid_smooth(spec1, stddev, inplace=False)
    assert np.allclose(spec1_smoothed, flux_smoothed)

    # Test in-place smoothing
    spec1_smoothed = trapezoid_smooth(spec1, stddev, inplace=True)
    assert np.allclose(spec1_smoothed, flux_smoothed)
