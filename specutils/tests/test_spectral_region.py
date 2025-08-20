import numpy as np
import astropy.units as u
from numpy.testing import assert_allclose
import pytest
from astropy.wcs import WCS

from specutils.spectra.spectrum import Spectrum
from specutils.spectra.spectral_region import SpectralRegion
from specutils.manipulation import extract_region

@pytest.fixture
def frequency_spectrum():
    # Create basic frequency WCS
    w = WCS(naxis=1)
    w.wcs.crval = [1.410e9]     # starting frequency (Hz)
    w.wcs.cdelt = [1.0e6]       # 1 MHz per channel
    w.wcs.crpix = [1]           # reference pixel
    w.wcs.cunit = ['Hz']
    w.wcs.restfrq = 1.420e9     # rest frequency in Hz

    # Build spectral axis and flux
    freqs = np.arange(1.410e9, 1.431e9, 1.0e6) * u.Hz
    flux = np.arange(1, len(freqs) + 1, dtype=float) * u.Jy

    return Spectrum(spectral_axis=freqs, flux=flux, wcs=w, velocity_convention='radio')


def test_extract_region_velocity_on_frequency_axis(frequency_spectrum):
    spec = frequency_spectrum

    # Define velocity range
    region = SpectralRegion(-500 * u.km / u.s, 500 * u.km / u.s)

    # Extract region with WCS preservation
    sub = extract_region(spec, region, preserve_wcs=True)

    # Determine expected frequency channels based on velocity condition
    velocities = spec.velocity.to(u.km / u.s)
    mask = (velocities >= -500 * u.km / u.s) & (velocities <= 500 * u.km / u.s)
    expected_freqs = spec.spectral_axis[mask]
    expected_flux = spec.flux[mask]

    # Assertions
    assert len(sub.spectral_axis) == len(expected_freqs)
    assert_allclose(sub.spectral_axis.to_value(u.Hz),
                    expected_freqs.to_value(u.Hz),
                    rtol=0, atol=1e-12)

    assert_allclose(sub.flux.to_value(u.Jy),
                    expected_flux.to_value(u.Jy),
                    rtol=0, atol=0)

    assert np.isclose(sub.wcs.wcs.crval[0], expected_freqs[0].value)
    assert np.isclose(sub.wcs.wcs.crpix[0], 1)
    assert np.isclose(sub.wcs.wcs.cdelt[0], spec.wcs.wcs.cdelt[0])
    assert sub.wcs.wcs.restfrq == spec.wcs.wcs.restfrq

def test_extract_region_drops_wcs_when_disabled(frequency_spectrum):
    spec = frequency_spectrum

    # Define velocity range
    region = SpectralRegion(-500 * u.km / u.s, 500 * u.km / u.s)

    # Extract region without WCS preservation
    sub = extract_region(spec, region, preserve_wcs=False)

    # Basic content check
    velocities = spec.velocity.to(u.km / u.s)
    mask = (velocities >= -500 * u.km / u.s) & (velocities <= 500 * u.km / u.s)
    expected_flux = spec.flux[mask]

    assert_allclose(sub.flux.to_value(u.Jy),
                    expected_flux.to_value(u.Jy),
                    rtol=0, atol=0)

    # Ensure WCS is removed
    assert not hasattr(sub, "wcs") or sub.wcs is None
