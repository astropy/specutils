"""Helpers for testing specutils objects in ASDF files.
These are similar to those in ``asdf_astropy.testing.helpers``.
"""
from astropy.tests.helper import assert_quantity_allclose
from numpy.testing import assert_allclose

__all__ = ["assert_spectrum1d_equal", "assert_spectrumlist_equal"]


def assert_spectrum1d_equal(a, b):
    """Equality test for use in ASDF unit tests for Spectrum1D."""
    __tracebackhide__ = True

    assert_quantity_allclose(a.flux, b.flux)
    assert_quantity_allclose(a.spectral_axis, b.spectral_axis)
    if a.uncertainty is None:
        assert b.uncertainty is None
    else:
        assert a.uncertainty.uncertainty_type == b.uncertainty.uncertainty_type
        assert_allclose(a.uncertainty.array, b.uncertainty.array)


def assert_spectrumlist_equal(a, b):
    """Equality test for use in ASDF unit tests for SpectrumList."""
    __tracebackhide__ = True

    assert len(a) == len(b)
    for x, y in zip(a, b):
        assert_spectrum1d_equal(x, y)
