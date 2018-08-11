import astropy.units as u
import numpy as np
import pytest

from ..spectra.spectrum1d import Spectrum1D
from ..spectra.spectrum_collection import SpectrumCollection


def test_create_spectrum_collection():
    spec = Spectrum1D(spectral_axis=np.linspace(0, 50, 50) * u.AA,
                      flux=np.random.randn(50) * u.Jy)
    spec1 = Spectrum1D(spectral_axis=np.linspace(0, 50, 25) * u.AA,
                       flux=np.random.randn(25) * u.Jy)

    sc = SpectrumCollection([spec, spec1], output_grid='coarse')

    assert len(sc) == 2
    assert isinstance(sc[0], Spectrum1D)
    assert sc.flux.shape == (2, 25)
    assert isinstance(sc.flux, u.Quantity)
    assert isinstance(sc.flux[0], u.Quantity)

    sc1 = SpectrumCollection([spec, spec1], output_grid='fine')

    with pytest.raises(Exception):
        SpectrumCollection([spec, spec1], output_grid='same')

    sc4 = SpectrumCollection([spec, spec1], output_grid=(0, 30, 1))
    sc5 = SpectrumCollection([spec, spec1], output_grid=[0., 1.42857143, 2.85714286, 4.28571429, 5.71428571, 7.14285714, 8.57142857, 10.])

    assert sc1.flux.shape == (2, 50)
    assert sc4.flux.shape == (2, 30)
    assert sc5.flux.shape == (2, 8)
