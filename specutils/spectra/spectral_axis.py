import warnings

import astropy.units as u
import numpy as np

from .spectral_coordinate import SpectralCoord


__all__ = ['SpectralAxis']

# We don't want to run doctests in the docstrings we inherit from Quantity
__doctest_skip__ = ['SpectralAxis.*']

class SpectralAxis(SpectralCoord):
    """
    """

    def __new__():
        pass

    @static_method
    def _edges_from_centers():
        pass

    @static_method
    def centers_from_edges():
        pass

    @property
    def bin_edges(self):
        pass

    @property
    def bin_centers(self):
        pass
