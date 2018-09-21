import pytest
import numpy as np

import astropy.units as u
from astropy.tests.helper import quantity_allclose

from ..spectra import Spectrum1D, SpectralRegion


from ..manipulation import noise_region_uncertainty


def pure_noise_spectrum(amplitude=1*u.mJy):
    np.random.seed(42)

    lamb = np.linspace(4000, 10000, 1000)*u.AA
    flux = np.random.randn(1000) * amplitude
    return Spectrum1D(spectral_axis=lamb, flux=flux)


def test_noise_region_unc():
    nspec = pure_noise_spectrum()

    region = SpectralRegion(5000*u.AA, 6000*u.AA)

    unc_spec = noise_region_uncertainty(nspec, region)

    assert quantity_allclose(unc_spec.uncertainty, 1*u.mJy, atol=.1) #this is a guess at the atol... need to finalize with real data


def test_noise_region_unc_multi():
    nspec = pure_noise_spectrum()

    multi_region = SpectralRegion([(5000*u.AA, 6000*u.AA), (900*u.pixel, 980*u.pixel)])

    unc_spec = noise_region_uncertainty(nspec, multi_region)

    assert quantity_allclose(unc_spec.uncertainty, 1*u.mJy, atol=.01) #this is a guess at the atol... need to finalize with real data
