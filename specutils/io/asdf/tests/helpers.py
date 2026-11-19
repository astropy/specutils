"""Helpers for testing specutils objects in ASDF files.
These are similar to those in ``asdf_astropy.testing.helpers``.
"""
import astropy.units as u
from asdf_astropy.testing.helpers import assert_frame_equal
from astropy.coordinates import SkyCoord
from astropy.tests.helper import assert_quantity_allclose
from numpy.testing import assert_allclose, assert_array_equal

__all__ = ["assert_spectral_axis_equal", "assert_spectrum_equal", "assert_spectrumlist_equal"]


def assert_spectral_axis_equal(a, b):
    """Equality test for use in ASDF unit tests for SpectralAxis."""
    __tracebackhide__ = True

    assert type(a) is type(b)
    assert_quantity_allclose(a.quantity, b.quantity)
    assert_frame_equal(a.target, b.target)
    # The observer may come back in a different representation; compare
    # its position and velocity instead
    assert (a.observer is None) == (b.observer is None)
    if a.observer is not None:
        a_obs, b_obs = SkyCoord(a.observer).icrs, SkyCoord(b.observer).icrs
        assert_quantity_allclose(a_obs.cartesian.xyz, b_obs.cartesian.xyz)
        assert_quantity_allclose(a_obs.velocity.d_xyz, b_obs.velocity.d_xyz)
    assert a.medium == b.medium
    assert a.frame == b.frame
    assert a.obstime == b.obstime
    if a.location is None:
        assert b.location is None
    else:
        assert_quantity_allclose(a.location.geocentric, b.location.geocentric)
    assert_quantity_allclose(a.radial_velocity, b.radial_velocity, atol=1e-6 * u.km / u.s)
    assert a.doppler_rest == b.doppler_rest
    assert a.doppler_convention == b.doppler_convention


def assert_spectrum_equal(a, b):
    """Equality test for use in ASDF unit tests for Spectrum."""
    __tracebackhide__ = True

    assert_quantity_allclose(a.flux, b.flux)
    assert_spectral_axis_equal(a.spectral_axis, b.spectral_axis)

    if a.uncertainty is None:
        assert b.uncertainty is None
    else:
        assert a.uncertainty.uncertainty_type == b.uncertainty.uncertainty_type
        assert_allclose(a.uncertainty.array, b.uncertainty.array)

    if a.mask is None:
        assert b.mask is None
    else:
        assert_array_equal(a.mask, b.mask)


def assert_spectrumlist_equal(a, b):
    """Equality test for use in ASDF unit tests for SpectrumList."""
    __tracebackhide__ = True

    assert len(a) == len(b)
    for x, y in zip(a, b):
        assert_spectrum_equal(x, y)
