import astropy.units as u
import astropy.wcs as fitswcs
import gwcs
import numpy as np
import pytest

from ..spectra.spectrum1d import Spectrum1D


def test_spectral_axes():
    flux1 = (np.random.sample(49) * 100).astype(int)
    flux2 = (np.random.sample(49) * 100).astype(int)

    flux3 = flux1 + flux2

    spec1 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                       flux=flux1)
    spec2 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                       flux=flux2)

    spec3 = spec1 + spec2

    assert np.allclose(spec3.flux.value, flux3)
