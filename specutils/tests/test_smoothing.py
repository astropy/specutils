import numpy as np
import pytest
from astropy import convolution
from scipy.signal import medfilt
import astropy.units as u
from astropy.nddata import StdDevUncertainty, VarianceUncertainty, InverseVariance

from ..spectra import Spectrum1D
from ..manipulation.smoothing import (convolution_smooth, box_smooth,
                                      gaussian_smooth, trapezoid_smooth,
                                      median_smooth)


def compare_flux(flux_smooth1, flux_smooth2, flux_original, rtol=0.01):
    """
    There are two things to compare for each set of smoothing:

    1. Compare the smoothed flux from the astropy machinery vs
       the smoothed flux from specutils.  This is done by
       comparing flux_smooth1 and flux_smooth2.

    2. Next we want to compare the smoothed flux to the original
       flux.  This is a little more difficult as smoothing will
       make a difference for median filter, but less so for
       convolution based smoothing if the kernel is normalized
       (area under the kernel = 1).

       In this second case the rtol (relative tolerance) is used
       judiciously.

    """

    # Compare, element by element, the two smoothed fluxes.
    # Workaround for astropy dev (4.2) which seems to return quantities in
    #  convolutions
    if isinstance(flux_smooth1, u.Quantity):
        flux_smooth1 = flux_smooth1.value
    if isinstance(flux_smooth2, u.Quantity):
        flux_smooth2 = flux_smooth2.value

    assert np.allclose(flux_smooth1, flux_smooth2)

    # Compare the total spectral flux of the smoothed to the original.
    assert np.allclose(sum(flux_smooth1), sum(flux_original), rtol=rtol)


def test_smooth_custom_kernel(simulated_spectra):
    """
    Test CustomKernel smoothing with correct parameters.
    """

    # Create the original spectrum
    spec1 = simulated_spectra.s1_um_mJy_e1
    flux_original = spec1.flux

    # Create a custom kernel (some weird asymmetric-ness)
    numpy_kernel = np.array([0.5, 1, 2, 0.5, 0.2])
    numpy_kernel = numpy_kernel / np.sum(numpy_kernel)

    custom_kernel = convolution.CustomKernel(numpy_kernel)
    flux_smoothed_astropy = convolution.convolve(flux_original, custom_kernel)

    # Calculate the custom smoothed
    spec1_smoothed = convolution_smooth(spec1, custom_kernel)
    compare_flux(spec1_smoothed.flux.value, flux_smoothed_astropy, flux_original.value)


@pytest.mark.parametrize("width", [1, 2.3])
def test_smooth_box_good(simulated_spectra, width):
    """
    Test Box1DKernel smoothing with correct parameters.

    Width values need to be a number greater than 0.
    """

    # Create the original spectrum
    spec1 = simulated_spectra.s1_um_mJy_e1
    flux_original = spec1.flux

    # Calculate the smoothed flux using Astropy
    box_kernel = convolution.Box1DKernel(width)
    flux_smoothed_astropy = convolution.convolve(flux_original, box_kernel)

    # Calculate the box smoothed
    spec1_smoothed = box_smooth(spec1, width)
    compare_flux(spec1_smoothed.flux.value, flux_smoothed_astropy, flux_original.value)

    # Check the input and output units
    assert spec1.wavelength.unit == spec1_smoothed.wavelength.unit
    assert spec1.flux.unit == spec1_smoothed.flux.unit
    assert len(spec1.meta) == len(spec1_smoothed.meta)


@pytest.mark.parametrize("width", [-1, 0, 'a'])
def test_smooth_box_bad(simulated_spectra, width):
    """
    Test Box1DKernel smoothing with incorrect parameters.

    Width values need to be a number greater than 0.
    """

    #  Create the spectrum
    spec1 = simulated_spectra.s1_um_mJy_e1

    # Test bad input parameters
    with pytest.raises(ValueError):
        box_smooth(spec1, width)


@pytest.mark.parametrize("stddev", [1, 2.3])
def test_smooth_gaussian_good(simulated_spectra, stddev):
    """
    Test Gaussian1DKernel smoothing with correct parameters.

    Standard deviation values need to be a number greater than 0.
    """

    #  Create the spectrum
    spec1 = simulated_spectra.s1_um_mJy_e1
    flux_original = spec1.flux

    # Calculate the smoothed flux using Astropy
    gaussian_kernel = convolution.Gaussian1DKernel(stddev)
    flux_smoothed_astropy = convolution.convolve(flux_original, gaussian_kernel)

    # Test gaussian smoothing
    spec1_smoothed = gaussian_smooth(spec1, stddev)
    compare_flux(spec1_smoothed.flux.value, flux_smoothed_astropy, flux_original.value, rtol=0.02)

    # Check the input and output units
    assert spec1.wavelength.unit == spec1_smoothed.wavelength.unit
    assert spec1.flux.unit == spec1_smoothed.flux.unit
    assert len(spec1.meta) == len(spec1_smoothed.meta)


@pytest.mark.parametrize("stddev", [-1, 0, 'a'])
def test_smooth_gaussian_bad(simulated_spectra, stddev):
    """
    Test MexicanHat1DKernel smoothing with incorrect parameters.

    Standard deviation values need to be a number greater than 0.
    """

    #  Create the spectrum
    spec1 = simulated_spectra.s1_um_mJy_e1

    # Test bad input paramters
    with pytest.raises(ValueError):
        gaussian_smooth(spec1, stddev)


@pytest.mark.parametrize("stddev", [1, 2.3])
def test_smooth_trapezoid_good(simulated_spectra, stddev):
    """
    Test Trapezoid1DKernel smoothing with correct parameters.

    Standard deviation values need to be a number greater than 0.
    """

    #  Create the spectrum
    spec1 = simulated_spectra.s1_um_mJy_e1
    flux_original = spec1.flux

    # Create the flux_smoothed which is what we want to compare to
    trapezoid_kernel = convolution.Trapezoid1DKernel(stddev)
    flux_smoothed_astropy = convolution.convolve(flux_original, trapezoid_kernel)

    # Test trapezoid smoothing
    spec1_smoothed = trapezoid_smooth(spec1, stddev)
    compare_flux(spec1_smoothed.flux.value, flux_smoothed_astropy, flux_original.value)

    # Check the input and output units
    assert spec1.wavelength.unit == spec1_smoothed.wavelength.unit
    assert spec1.flux.unit == spec1_smoothed.flux.unit
    assert len(spec1.meta) == len(spec1_smoothed.meta)


@pytest.mark.parametrize("stddev", [-1, 0, 'a'])
def test_smooth_trapezoid_bad(simulated_spectra, stddev):
    """
    Test Trapezoid1DKernel smoothing with incorrect parameters.

    Standard deviation values need to be a number greater than 0.
    """

    #  Create the spectrum
    spec1 = simulated_spectra.s1_um_mJy_e1

    # Test bad parameters
    with pytest.raises(ValueError):
        trapezoid_smooth(spec1, stddev)


@pytest.mark.parametrize("width", [1, 3, 9])
def test_smooth_median_good(simulated_spectra, width):
    """
    Test Median smoothing with correct parameters.

    Width values need to be a number greater than 0.
    """

    #  Create the spectrum
    spec1 = simulated_spectra.s1_um_mJy_e1
    flux_original = spec1.flux

    # Create the flux_smoothed which is what we want to compare to
    flux_smoothed_astropy = medfilt(flux_original, width)

    # Test median smoothing
    spec1_smoothed = median_smooth(spec1, width)
    compare_flux(spec1_smoothed.flux.value, flux_smoothed_astropy, flux_original.value, rtol=0.15)

    # Check the input and output units
    assert spec1.wavelength.unit == spec1_smoothed.wavelength.unit
    assert spec1.flux.unit == spec1_smoothed.flux.unit
    assert len(spec1.meta) == len(spec1_smoothed.meta)


@pytest.mark.parametrize("width", [-1, 0, 'a'])
def test_smooth_median_bad(simulated_spectra, width):
    """
    Test Median smoothing with incorrect parameters.

    Width values need to be a number greater than 0.
    """

    #  Create the spectrum
    spec1 = simulated_spectra.s1_um_mJy_e1

    # Test bad parameters
    with pytest.raises(ValueError):
        median_smooth(spec1, width)


def test_smooth_custom_kernel_uncertainty(simulated_spectra):
    """
    Test CustomKernel smoothing with correct parameters.
    """

    np.random.seed(42)

    # Create a custom kernel (some weird asymmetric-ness)
    numpy_kernel = np.array([0.5, 1, 2, 0.5, 0.2])
    numpy_kernel = numpy_kernel / np.sum(numpy_kernel)
    custom_kernel = convolution.CustomKernel(numpy_kernel)

    spec1 = simulated_spectra.s1_um_mJy_e1
    uncertainty = np.abs(np.random.random(spec1.flux.shape))

    # Test StdDevUncertainty
    spec1.uncertainty = StdDevUncertainty(uncertainty)

    spec1_smoothed = convolution_smooth(spec1, custom_kernel)
    tt = convolution.convolve(1/(spec1.uncertainty.array**2), custom_kernel)
    uncertainty_smoothed_astropy = 1/np.sqrt(tt)

    assert np.allclose(spec1_smoothed.uncertainty.array, uncertainty_smoothed_astropy)

    # Test VarianceUncertainty
    spec1.uncertainty = VarianceUncertainty(uncertainty)

    spec1_smoothed = convolution_smooth(spec1, custom_kernel)
    uncertainty_smoothed_astropy = 1/convolution.convolve(1/spec1.uncertainty.array, custom_kernel)
    assert np.allclose(spec1_smoothed.uncertainty.array, uncertainty_smoothed_astropy)

    # Test InverseVariance
    spec1.uncertainty = InverseVariance(uncertainty)

    spec1_smoothed = convolution_smooth(spec1, custom_kernel)
    uncertainty_smoothed_astropy = convolution.convolve(spec1.uncertainty.array, custom_kernel)
    assert np.allclose(spec1_smoothed.uncertainty.array, uncertainty_smoothed_astropy)


def test_medfilt_dtype():
    """
    Check that dtypes not supported by scipy's medfilt method are
    converted to supported dtypes before they are passed to scipy.
    """
    flux = [1, 2, 3] * u.Jy

    for dtype in ('>f4', np.float32, 'i4', int):
        spectrum = Spectrum1D(flux=flux.astype(dtype), spectral_axis=[1, 2, 3] * u.um)
        median_smooth(spectrum, 1)
