import astropy.units as u
import numpy as np
from numpy.testing import assert_allclose

from ..spectra.spectrum1d import Spectrum1D


def test_spectral_axes():
    spec1 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                       flux=np.random.sample(49) * 100)

    sliced_spec1 = spec1[0]

    assert isinstance(sliced_spec1, Spectrum1D)
    #assert sliced_spec1.wcs == spec1.wcs
    assert_allclose(sliced_spec1.wcs.pixel_to_world(np.arange(10)), spec1.wcs.pixel_to_world(np.arange(10)))

    flux2 = np.random.sample((10, 49)) * 100

    spec2 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                       flux=flux2)

    sliced_spec2 = spec2[0]

    assert isinstance(sliced_spec2, Spectrum1D)
    #assert sliced_spec2.wcs == spec2.wcs
    assert_allclose(sliced_spec2.wcs.pixel_to_world(np.arange(10)), spec2.wcs.pixel_to_world(np.arange(10)))
    assert sliced_spec2.flux.shape[0] == 49


def test_slicing():

    # Create the initial spectrum
    spec = Spectrum1D(spectral_axis=np.arange(10) * u.um, flux=2*np.arange(10)*u.Jy)

    # Slice it.
    sub_spec = spec[4:8]

    # Check basic spectral_axis property
    assert sub_spec.spectral_axis.unit == u.um
    assert np.allclose(sub_spec.spectral_axis.value, np.array([4, 5, 6, 7]))
    assert np.allclose(sub_spec.flux.value, np.array([8, 10, 12, 14]))

    assert sub_spec.wavelength.unit == u.AA
    assert np.allclose(sub_spec.wavelength.value, np.array([40000., 50000., 60000., 70000.]))

    assert sub_spec.frequency.unit == u.GHz
    assert np.allclose(sub_spec.frequency.value, np.array([74948.1145, 59958.4916, 49965.40966667, 42827.494]))

    # Do it a second time to confirm the original was not modified.
    sub_spec2 = spec[1:5]

    # Check basic spectral_axis property
    assert sub_spec2.spectral_axis.unit == u.um
    assert np.allclose(sub_spec2.spectral_axis.value, np.array([1, 2, 3, 4]))
    assert np.allclose(sub_spec2.flux.value, np.array([2, 4, 6, 8]))

    assert sub_spec2.wavelength.unit == u.AA
    assert np.allclose(sub_spec2.wavelength.value, np.array([10000., 20000., 30000., 40000.]))

    assert sub_spec2.frequency.unit == u.GHz
    assert np.allclose(sub_spec2.frequency.value, np.array([299792.458, 149896.229,  99930.81933333,  74948.1145]))

    # Going to repeat these to make sure the original spectrum was
    # not modified in some way
    assert spec.spectral_axis.unit == u.um
    assert np.allclose(spec.spectral_axis.value, np.array(np.arange(10)))
    assert np.allclose(spec.flux.value, np.array(2*np.arange(10)))

    assert spec.wavelength.unit == u.AA
    assert np.allclose(spec.wavelength.value, np.array(10000*np.arange(10)))

    assert sub_spec.frequency.unit == u.GHz
    assert np.allclose(sub_spec.frequency.value, np.array([74948.1145, 59958.4916, 49965.40966667, 42827.494]))

