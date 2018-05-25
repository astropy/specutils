import astropy.units as u
import numpy as np

from .spectral_examples import spectral_examples
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


def test_add_spectra(spectral_examples):

    # Get the numpy array of data
    flux1 = define_spectra.s1_um_mJy_e1_flux
    flux2 = define_spectra.s1_um_mJy_e2_flux
    flux3 = flux1 + flux2

    # Calculate using the spectrum1d/nddata code
    spec3 = define_spectra.s1_um_mJy_e1 + define_spectra.s1_um_mJy_e2

    assert np.allclose(spec3.flux.value, flux3)
