import astropy.units as u
import astropy.wcs as fitswcs
from astropy.tests.helper import quantity_allclose
import numpy as np
from numpy.testing import assert_allclose

from ..spectra.spectrum1d import Spectrum1D


def test_spectral_axes():
    spec1 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                       flux=np.random.sample(49) * 100 * u.Jy)

    sliced_spec1 = spec1[0:2]

    assert isinstance(sliced_spec1, Spectrum1D)
    assert_allclose(sliced_spec1.wcs.pixel_to_world(0), spec1.wcs.pixel_to_world(0))

    flux2 = np.random.sample((10, 49)) * 100

    spec2 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                       flux=flux2 * u.Jy)

    sliced_spec2 = spec2[0]

    assert isinstance(sliced_spec2, Spectrum1D)
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
    assert np.allclose(sub_spec2.frequency.value, np.array([299792.458, 149896.229, 99930.81933333, 74948.1145]))

    # Going to repeat these to make sure the original spectrum was
    # not modified in some way
    assert spec.spectral_axis.unit == u.um
    assert np.allclose(spec.spectral_axis.value, np.array(np.arange(10)))
    assert np.allclose(spec.flux.value, np.array(2*np.arange(10)))

    assert spec.wavelength.unit == u.AA
    assert np.allclose(spec.wavelength.value, np.array(10000*np.arange(10)))

    assert sub_spec.frequency.unit == u.GHz
    assert np.allclose(sub_spec.frequency.value, np.array([74948.1145, 59958.4916, 49965.40966667, 42827.494]))


def test_slicing_with_fits():
    my_wcs = fitswcs.WCS(header={'CDELT1': 1, 'CRVAL1': 6562.8, 'CUNIT1': 'Angstrom',
                                 'CTYPE1': 'WAVE', 'RESTFRQ': 1400000000, 'CRPIX1': 25})

    spec = Spectrum1D(flux=[5, 6, 7, 8, 9, 10] * u.Jy, wcs=my_wcs)
    spec_slice = spec[1:5]

    assert isinstance(spec_slice, Spectrum1D)
    assert spec_slice.flux.size == 4
    assert quantity_allclose(spec_slice.wcs.pixel_to_world([0, 1, 2, 3]),
                             spec.wcs.pixel_to_world([1, 2, 3, 4]))


def test_slicing_multidim():
    spec = Spectrum1D(spectral_axis=np.arange(10) * u.AA,
                      flux=np.random.sample((5, 10)) * u.Jy,
                      mask=np.random.sample((5, 10)) > 0.5)

    spec1 = spec[0]
    spec2 = spec[1:3]
    spec3 = spec[..., 4:7]

    assert spec1.flux[0] == spec.flux[0][0]
    assert quantity_allclose(spec1.spectral_axis, spec.spectral_axis)
    assert spec.flux.shape[1:] == spec1.flux.shape

    assert quantity_allclose(spec2.flux, spec.flux[1:3])
    assert quantity_allclose(spec2.spectral_axis, spec.spectral_axis)

    assert spec1.mask[0] == spec.mask[0][0]
    assert spec1.mask.shape == (10,)

    assert quantity_allclose(spec3.spectral_axis, spec.spectral_axis[4:7])
    assert quantity_allclose(spec3.wcs.pixel_to_world([0,1,2]), spec3.spectral_axis[0:3])
