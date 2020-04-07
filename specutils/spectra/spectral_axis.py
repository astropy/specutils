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

    def __new__(cls, value, unit=None, observer=None, target=None,
                 radial_velocity=None, redshift=None, doppler_rest=None,
                 doppler_convention=None, bin_specification="centers",
                 **kwargs):

        # Convert to bin centers if bin edges were given, since SpectralCoord
        # only accepts centers
        if bin_specification == "edges":
            value = SpectralAxis._centers_from_edges(value)

        obj = super().__new__(cls, value, unit=unit, observer=observer,
                              target=target, radial_velocity=radial_velocity,
                              redshift=redshift, doppler_rest=doppler_rest,
                              doppler_convention=doppler_convention, **kwargs)

        #obj.bin_edges = cls._edges_from_centers(value)

        return obj

    def __quantity_subclass__(self, unit):
        """:wq

        Overridden by subclasses to change what kind of view is
        created based on the output unit of an operation.
        """
        return SpectralAxis, True

    @staticmethod
    def _edges_from_centers(centers):
        a = np.insert(centers, 0, 2*centers[0]-centers[1])
        b = np.append(centers, 2*centers[-1]-centers[-2])
        edges = (a + b) / 2
        return edges

    @staticmethod
    def _centers_from_edges(edges):
        return (edges[1:] + edges[:-1]) / 2

    @property
    def bin_edges(self):
        return self._edges_from_centers(self.value)
