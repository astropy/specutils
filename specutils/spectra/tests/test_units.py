from astropy import units as u
import numpy as np
import pytest

from ..spectrum1d import Spectrum1D


def test_flux_unit_conversion():
    # By default the flux units should be set to Jy
    s = Spectrum1D(flux=np.array([26.0, 44.5]), spectral_axis=np.array([400, 500]) * u.nm)
    assert np.all(s.flux == np.array([26.0, 44.5]) * u.Jy)
    assert s.flux.unit == u.Jy

    # Simple Unit Conversion
    s = Spectrum1D(flux=np.array([26.0, 44.5]) * u.Jy, spectral_axis=np.array([400, 500])*u.nm)
    converted_value = s.to_flux(unit=u.uJy)[0]
    assert ((26.0 * u.Jy).to(u.uJy) == converted_value)

    # Make sure incompatible units raise UnitConversionError
    with pytest.raises(u.UnitConversionError):
        converted_value = s.to_flux(unit=u.m)

    # Pass custom equivalencies
    s = Spectrum1D(flux=np.array([26.0, 44.5]) * u.Jy, spectral_axis=np.array([400, 500]) * u.nm)
    eq = [[u.Jy, u.m,
          lambda x: np.full_like(np.array(x), 1000.0, dtype=np.double),
          lambda x: np.full_like(np.array(x), 0.001, dtype=np.double)]]
    converted_value = s.to_flux(unit=u.m, equivalencies=eq)[0]
    assert 1000.0 * u.m == converted_value

    # Check if suppressing the unit conversion works
    s = Spectrum1D(flux=np.array([26.0, 44.5]) * u.Jy, spectral_axis=np.array([400, 500]) * u.nm)
    s.to_flux("uJy", suppress_conversion=True)
    assert s.flux[0] == 26.0 * u.uJy


def test_spectral_axis_unit_conversion():
    # By default the spectral axis units should be set to angstroms
    s = Spectrum1D(flux=np.array([26.0, 44.5]) * u.Jy, spectral_axis=np.array([400, 500]))
    assert np.all(s.spectral_axis == np.array([400, 500]) * u.angstrom)
    assert s.spectral_axis.unit == u.angstrom

    # Simple Unit Conversion
    s = Spectrum1D(flux=np.array([26.0, 44.5]) * u.Jy, spectral_axis=np.array([400, 500]) * u.nm)

    new_s = s.with_spectral_unit("m")
    assert np.all(s.spectral_axis.to("m") == new_s.spectral_axis)
