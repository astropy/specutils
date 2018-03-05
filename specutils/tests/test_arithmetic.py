import astropy.units as u
import astropy.wcs as fitswcs
import gwcs
import numpy as np
import pytest

from ..spectra.spectrum1d import Spectrum1D


def test_spectral_axes():
    spec1 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                       flux=np.random.randn(49))
    spec2 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                       flux=np.random.randn(49))

    spec3 = spec1 + spec2

    # Arithmetic should defautly fail when the axes are different
    spec1 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                       flux=np.random.randn(49))
    spec2 = Spectrum1D(spectral_axis=np.arange(1100, 1150) * u.Angstrom,
                       flux=np.random.randn(49))
