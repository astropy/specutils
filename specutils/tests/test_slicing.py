import astropy.units as u
import numpy as np
from numpy.testing import assert_allclose

from ..spectra.spectrum1d import Spectrum1D


def test_spectral_axes():
    spec1 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                       flux=np.random.sample(49) * 100)

    sliced_spec1 = spec1[0]

    assert isinstance(sliced_spec1, Spectrum1D)
    assert sliced_spec1.wcs == spec1.wcs
    assert_allclose(sliced_spec1.wcs.pixel_to_world(np.arange(10)), spec1.wcs.pixel_to_world(np.arange(10)))

    flux2 = np.random.sample((10, 49)) * 100

    spec2 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                       flux=flux2)

    sliced_spec2 = spec2[0]

    assert isinstance(sliced_spec2, Spectrum1D)
    assert sliced_spec2.wcs == spec2.wcs
    assert_allclose(sliced_spec2.wcs.pixel_to_world(np.arange(10)), spec2.wcs.pixel_to_world(np.arange(10)))
    assert sliced_spec2.flux.shape[0] == 49


def test_slicing():

    # Create the initial spectrum
    spec = Spectrum1D(spectral_axis=np.arange(10) * u.AA, flux=2*np.arange(10)*u.Jy)

    # Slice it.
    sub_spec = spec[4:8]

    assert sub_spec.spectral_axis.unit == u.AA
    assert np.allclose(sub_spec.spectral_axis.value, np.array([4, 5, 6, 7]))
    assert np.allclose(sub_spec.flux.value, np.array([8, 10, 12, 14]))
