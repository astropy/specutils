import astropy.units as u
import numpy as np
import pytest

from ..spectra.spectrum1d import Spectrum1D
from ..spectra.spectrum_collection import SpectrumCollection


def test_create_spectrum_collection():
    spec = Spectrum1D(spectral_axis=np.linspace(0, 50, 50) * u.AA,
                      flux=np.random.randn(50) * u.Jy)
    spec1 = Spectrum1D(spectral_axis=np.linspace(20, 60, 50) * u.AA,
                       flux=np.random.randn(50) * u.Jy)

    sc = SpectrumCollection([spec, spec1])

    assert len(sc) == 2
    assert isinstance(sc[0], Spectrum1D)
    assert isinstance(sc.flux, u.Quantity)
    assert sc.flux.shape == (2, 50)
    assert isinstance(sc.flux[0], u.Quantity)
    assert sc.flux[0].size == 1

    # Test to ensure that the shape requirement is enforced
    with pytest.raises(Exception):
        spec = Spectrum1D(spectral_axis=np.linspace(0, 50, 25) * u.AA,
                          flux=np.random.randn(25) * u.Jy)
        spec1 = Spectrum1D(spectral_axis=np.linspace(20, 60, 50) * u.AA,
                           flux=np.random.randn(50) * u.Jy)

        sc = SpectrumCollection([spec, spec1])

