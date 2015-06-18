from specutils import Spectrum1D
from astropy import units as u
import numpy as np
from numpy import testing as nptesting


def test_spectrum1d_fromarray_quantity():
    test_spec = Spectrum1D.from_array(np.arange(3000, 9000) * u.angstrom,
                          np.random.random(6000))
    wavelength = test_spec.wavelength
    assert hasattr(test_spec, 'wavelength')
    assert test_spec.dispersion_unit == u.angstrom

def test_spectrum1d_fromarray_quantity2():
    test_spec = Spectrum1D.from_array(np.arange(3000, 9000) * u.angstrom,
                          np.random.random(6000), dispersion_unit='nm')

    assert hasattr(test_spec, 'wavelength')

    assert test_spec.dispersion_unit == u.nm
    nptesting.assert_allclose(test_spec.wavelength.value,
                               np.arange(3000, 9000) / 10.)

def test_spectrum1d_flux1():
        test_spec = Spectrum1D.from_array(np.arange(3000, 9000) * u.angstrom,
                          np.random.random(6000) * u.erg/u.s,
                          )

        assert not hasattr(test_spec.data, 'unit')
        assert test_spec.flux.unit == u.erg / u.s
        assert test_spec.unit == u.erg / u.s

def test_spectrum1d_flux2():
        test_spec = Spectrum1D.from_array(np.arange(3000, 9000) * u.angstrom,
                          np.random.random(6000) * u.erg/u.s,
                          )

        assert not hasattr(test_spec.data, 'unit')
        assert test_spec.flux.unit == u.erg / u.s
        assert test_spec.unit == u.erg / u.s
        new_flux = np.random.random(6000) * u.W

        test_spec.flux = new_flux

        assert test_spec.flux.unit == u.erg / u.s

        nptesting.assert_allclose(new_flux.to(u.erg / u.s).value,
                                  test_spec.data)

        test_spec.flux[1000] = 1000 * u.erg/u.s

        nptesting.assert_allclose(test_spec.flux[1000], 1000 * u.erg/u.s)
