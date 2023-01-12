import astropy.units as u
from astropy.utils.decorators import lazyproperty
from astropy.coordinates import SpectralCoord
import numpy as np

__all__ = ['SpectralAxis']

# We don't want to run doctests in the docstrings we inherit from Quantity
__doctest_skip__ = ['SpectralAxis.*']


class SpectralAxis(SpectralCoord):
    """
    Coordinate object representing spectral values corresponding to a specific
    spectrum. Overloads SpectralCoord with additional information (currently
    only bin edges).

    Parameters
    ----------
    bin_specification: str, optional
        Must be "edges" or "centers". Determines whether specified axis values
        are interpreted as bin edges or bin centers. Defaults to "centers".
    """

    _equivalent_unit = SpectralCoord._equivalent_unit + (u.pixel,)

    def __new__(cls, value, *args, bin_specification="centers", **kwargs):

        # Enforce pixel axes are ascending
        if ((type(value) is u.quantity.Quantity) and
                (value.size > 1) and
                (value.unit is u.pix) and
                (value[-1] <= value[0])):
            raise ValueError("u.pix spectral axes should always be ascending")

        # Convert to bin centers if bin edges were given, since SpectralCoord
        # only accepts centers
        if bin_specification == "edges":
            bin_edges = value
            value = SpectralAxis._centers_from_edges(value)

        obj = super().__new__(cls, value, *args, **kwargs)

        if bin_specification == "edges":
            obj._bin_edges = bin_edges

        return obj

    @staticmethod
    def _edges_from_centers(centers, unit):
        """
        Calculates interior bin edges based on the average of each pair of
        centers, with the two outer edges based on extrapolated centers added
        to the beginning and end of the spectral axis.
        """
        a = np.insert(centers, 0, 2*centers[0] - centers[1])
        b = np.append(centers, 2*centers[-1] - centers[-2])
        edges = (a + b) / 2
        return edges*unit

    @staticmethod
    def _centers_from_edges(edges):
        """
        Calculates the bin centers as the average of each pair of edges
        """
        return (edges[1:] + edges[:-1]) / 2

    @lazyproperty
    def bin_edges(self):
        """
        Calculates bin edges if the spectral axis was created with centers
        specified.
        """
        if hasattr(self, '_bin_edges'):
            return self._bin_edges
        else:
            return self._edges_from_centers(self.value, self.unit)

    def with_observer_stationary_relative_to(self, frame,
                                             velocity=None,
                                             preserve_observer_frame=False):
        if self.unit is u.pixel:
            raise u.UnitsError("Cannot transform spectral coordinates in pixel units")
        return super().with_observer_stationary_relative_to(frame,
                                                            velocity=velocity,
                                                            preserve_observer_frame=preserve_observer_frame)

    def with_radial_velocity_shift(self, target_shift=None, observer_shift=None):
        if self.unit is u.pixel:
            raise u.UnitsError("Cannot transform spectral coordinates in pixel units")
        return super().with_radial_velocity_shift(target_shift=target_shift,
                                                  observer_shift=observer_shift)
