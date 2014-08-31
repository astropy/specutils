from specutils import Spectrum1D
from astropy import units as u
import numpy as np
from numpy import testing as nptesting


def test_spectrum1d_fromarray_quantity():
    test_spec = Spectrum1D.from_array(np.arange(3000, 9000) * u.angstrom,
                          np.random.random(6000))

    assert hasattr(test_spec, 'wavelength')

    assert test_spec.dispersion_unit == u.angstrom

def test_spectrum1d_fromarray_quantity2():
    test_spec = Spectrum1D.from_array(np.arange(3000, 9000) * u.angstrom,
                          np.random.random(6000), dispersion_unit='nm')

    assert hasattr(test_spec, 'wavelength')

    assert test_spec.dispersion_unit == u.nm
    nptesting.assert_allclose(test_spec.wavelength.value,
                               np.arange(3000, 9000) / 10.)



