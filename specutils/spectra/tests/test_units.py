from astropy import units as u
import numpy as np

from ..spectrum1d import Spectrum1D


def test_flux_unit_conversion():
    # Simple Unit Conversion
    s = Spectrum1D(flux=np.array([26.0]) * u.Jy, spectral_axis=np.array([400])*u.nm)
    converted_value = s.to_flux(unit=u.uJy)[0]
    assert ((26.0 * u.Jy).to(u.uJy) == converted_value)

    # # Try to convert a unitless flux
    # s = Spectrum1D(flux=np.array([26.0]), spectral_axis=np.array([400]) * u.nm)
    # converted_value = s.to_flux(unit=u.uJy)[0]
    # print(converted_value)
    # assert ((26.0 * u.Jy).to(u.uJy) == converted_value)



'''
def test_wave_unit_conversion():
    s = Spectrum1D(flux=np.array([26.0]) * u.Jy, spectral_axis=np.array([400]) * u.nm)
    converted_value = s.to_flux(unit=u.uJy)[0]
'''