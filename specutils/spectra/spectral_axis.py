import warnings

import astropy.units as u
from astropy.utils.decorators import lazyproperty
from astropy.coordinates import (SpectralCoord, EarthLocation, CartesianDifferential,
                                 frame_transform_graph)
from astropy.coordinates.spectral_coordinate import NoVelocityWarning, NoDistanceWarning
from astropy.time import Time
import numpy as np

from .spectral_frame import SpectralMedium, SPECTRAL_FRAMES, normalize_frame

__all__ = ['SpectralAxis']

# We don't want to run doctests in the docstrings we inherit from Quantity
__doctest_skip__ = ['SpectralAxis.*']

_ZERO_VELOCITY = CartesianDifferential([0, 0, 0] * u.km / u.s)


def _validate_obstime(obstime):
    if obstime is None:
        return None
    if not isinstance(obstime, Time):
        obstime = Time(obstime)
    return obstime


def _validate_location(location):
    if location is None:
        return None
    if not isinstance(location, EarthLocation):
        raise TypeError("location must be an astropy EarthLocation, "
                        f"got {type(location).__name__}")
    return location


def observer_for_frame(frame, location, obstime, target=None):
    """
    Construct the observer, at ``location`` and ``obstime``, whose rest frame
    is the spectral reference frame ``frame``.

    Parameters
    ----------
    frame : str
        A FITS ``SPECSYS`` code (see `~specutils.spectra.spectral_frame.SPECTRAL_FRAMES`).
    location : `~astropy.coordinates.EarthLocation`
        Where the spectrum was recorded.
    obstime : `~astropy.time.Time`
        When the spectrum was recorded.
    target : `~astropy.coordinates.SkyCoord`, optional
        The source, needed only for the ``'SOURCE'`` frame.

    Returns
    -------
    `~astropy.coordinates.BaseCoordinateFrame` or None
        A frame instance with position and velocity, or `None` if no
        observer can be constructed for ``frame`` (``'LOCALGRP'``,
        ``'CMBDIPOL'``, or ``'SOURCE'`` without a target).
    """
    frame = normalize_frame(frame)
    topocentric = location.get_gcrs(obstime)
    if frame == 'TOPOCENT':
        return topocentric

    if frame == 'SOURCE':
        if target is None:
            return None
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', NoVelocityWarning)
            probe = SpectralCoord(1 * u.AA, observer=topocentric, target=target)
            return probe.with_observer_stationary_relative_to(probe.target).observer

    frame_name = SPECTRAL_FRAMES[frame]
    if frame_name is None:
        return None

    frame_cls = frame_transform_graph.lookup_name(frame_name)
    frame_kwargs = {}
    if 'obstime' in frame_cls.frame_attributes:
        frame_kwargs['obstime'] = obstime
    stationary = topocentric.transform_to(frame_cls(**frame_kwargs))
    return stationary.realize_frame(
        stationary.cartesian.without_differentials().with_differentials(_ZERO_VELOCITY))


class SpectralAxis(SpectralCoord):
    """
    Coordinate object representing spectral values corresponding to a specific
    spectrum. Overloads SpectralCoord with additional information: bin edges,
    the medium the wavelengths are expressed in, the velocity reference frame
    of the values, and where and when the spectrum was recorded.

    Parameters
    ----------
    bin_specification: str, optional
        Must be "edges" or "centers". Determines whether specified axis values
        are interpreted as bin edges or bin centers. Defaults to "centers".
    medium : `~specutils.spectra.spectral_frame.SpectralMedium`, {'vacuum', 'air'} or dict, optional
        The medium the wavelengths are expressed in. An air medium may carry
        the refraction formula and the air conditions it assumes.
    frame : str, optional
        FITS ``SPECSYS`` code for the reference frame the spectral values are
        measured in, e.g. ``'TOPOCENT'``, ``'BARYCENT'`` or ``'SOURCE'`` (the
        rest frame of the source). See
        `~specutils.spectra.spectral_frame.SPECTRAL_FRAMES`. The
        ``radial_velocity`` and ``redshift`` of the axis are the velocity of
        the source relative to an observer at rest in this frame.
    obstime : `~astropy.time.Time` or str, optional
        The mid-point of the observation.
    location : `~astropy.coordinates.EarthLocation`, optional
        Where the spectrum was recorded. If ``frame``, ``obstime`` and
        ``location`` are all given and no ``observer`` is, an observer at that
        location and time, at rest in ``frame``, is constructed.
    """

    _equivalent_unit = SpectralCoord._equivalent_unit + (u.pixel,)

    _metadata_attributes = ('medium', 'frame', 'obstime', 'location')

    def __new__(cls, value, *args, bin_specification="centers", medium=None,
                frame=None, obstime=None, location=None, **kwargs):

        # Enforce pixel axes are ascending
        if ((type(value) is u.quantity.Quantity) and
                (value.size > 1) and
                (value.unit is u.pix) and
                (value[-1] <= value[0])):
            raise ValueError("u.pix spectral axes should always be ascending")

        medium = SpectralMedium.from_input(medium)
        frame = normalize_frame(frame)
        obstime = _validate_obstime(obstime)
        location = _validate_location(location)

        # Inherit metadata from an existing SpectralAxis unless overridden
        if medium is None:
            medium = getattr(value, 'medium', None)
        if frame is None:
            frame = getattr(value, 'frame', None)
        if obstime is None:
            obstime = getattr(value, 'obstime', None)
        if location is None:
            location = getattr(value, 'location', None)

        observer = kwargs.get('observer')
        if observer is None:
            observer = getattr(value, 'observer', None)
        target = kwargs.get('target')
        if target is None:
            target = getattr(value, 'target', None)

        if (observer is None and frame is not None and location is not None
                and obstime is not None):
            observer = observer_for_frame(frame, location, obstime, target=target)
            if observer is not None:
                kwargs['observer'] = observer

        if (observer is not None and target is not None and
                (kwargs.get('radial_velocity') is not None or
                 kwargs.get('redshift') is not None)):
            raise ValueError("Cannot specify radial_velocity or redshift when both an "
                             "observer and a target are defined (the observer may have "
                             "been constructed from frame, location and obstime). Set "
                             "the velocity of the source on the target instead.")

        # Convert to bin centers if bin edges were given, since SpectralCoord
        # only accepts centers
        if bin_specification == "edges":
            bin_edges = value
            value = SpectralAxis._centers_from_edges(value)

        # A target or observer without a distance is taken to be very distant
        # (a source) or in the solar system (an observer), which is what we
        # want for positions read from headers; no need to warn about it.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', NoDistanceWarning)
            obj = super().__new__(cls, value, *args, **kwargs)

        if bin_specification == "edges":
            obj._bin_edges = bin_edges
        elif isinstance(value, SpectralAxis) and hasattr(value, '_bin_edges'):
            obj._bin_edges = value._bin_edges

        obj._medium = medium
        obj._frame = frame
        obj._obstime = obstime
        obj._location = location
        obj._check_medium_unit()

        return obj

    def __array_finalize__(self, obj):
        super().__array_finalize__(obj)
        self._medium = getattr(obj, '_medium', None)
        self._frame = getattr(obj, '_frame', None)
        self._obstime = getattr(obj, '_obstime', None)
        self._location = getattr(obj, '_location', None)

    def _check_medium_unit(self):
        if (self._medium is not None and self._medium.is_air and self.unit is not u.pixel
                and not self.unit.is_equivalent(u.m)):
            raise u.UnitsError("An air medium can only be set on a wavelength axis, "
                               f"not one in units of {self.unit}")

    @property
    def medium(self):
        """
        The `~specutils.spectra.spectral_frame.SpectralMedium` the wavelengths
        are expressed in, or `None` if unknown.
        """
        return self._medium

    @property
    def frame(self):
        """
        FITS ``SPECSYS`` code of the reference frame the spectral values are
        measured in, or `None` if unknown.
        """
        return self._frame

    @property
    def obstime(self):
        """The mid-point of the observation as a `~astropy.time.Time`, or `None`."""
        return self._obstime

    @property
    def location(self):
        """Where the spectrum was recorded, as an `~astropy.coordinates.EarthLocation`."""
        return self._location

    @property
    def _metadata(self):
        """Dictionary of the metadata attributes this axis carries."""
        return {name: getattr(self, name) for name in self._metadata_attributes}

    def replicate(self, value=None, unit=None, observer=None, target=None,
                  radial_velocity=None, redshift=None, doppler_convention=None,
                  doppler_rest=None, copy=False, medium=None, frame=None,
                  obstime=None, location=None):
        """
        Return a replica of the `SpectralAxis`, optionally changing the
        values or attributes. See `~astropy.coordinates.SpectralCoord.replicate`
        for the meaning of the arguments; ``medium``, ``frame``, ``obstime``
        and ``location`` are carried over from this axis unless given.
        """
        new = super().replicate(value=value, unit=unit, observer=observer,
                                target=target, radial_velocity=radial_velocity,
                                redshift=redshift, doppler_convention=doppler_convention,
                                doppler_rest=doppler_rest, copy=copy)
        # SpectralCoord.replicate turns an unset radial velocity into an
        # explicit zero; keep it unset so that the replica behaves the same.
        if (radial_velocity is None and redshift is None and self._radial_velocity is None
                and (new.observer is None or new.target is None)):
            new._radial_velocity = None
        new._medium = (SpectralMedium.from_input(medium) if medium is not None
                       else self._medium)
        new._frame = normalize_frame(frame) if frame is not None else self._frame
        new._obstime = _validate_obstime(obstime) if obstime is not None else self._obstime
        new._location = (_validate_location(location) if location is not None
                         else self._location)
        new._check_medium_unit()
        return new

    def to(self, unit, equivalencies=[], doppler_rest=None, doppler_convention=None):
        """
        Return a new `SpectralAxis` with the specified unit. See
        `~astropy.coordinates.SpectralQuantity.to`.

        Air wavelengths can only be converted to other wavelength units, since
        the relation between wavelength and frequency, energy or velocity only
        holds in vacuum. Convert to vacuum first.
        """
        unit = u.Unit(unit)
        if (self._medium is not None and self._medium.is_air
                and self.unit is not u.pixel and not unit.is_equivalent(u.m)):
            raise u.UnitConversionError(
                f"Cannot convert air wavelengths to '{unit}': the relation between "
                f"wavelength and {unit.physical_type} only holds in vacuum. Convert "
                "the spectrum to vacuum wavelengths first.")
        return super().to(unit, equivalencies=equivalencies, doppler_rest=doppler_rest,
                          doppler_convention=doppler_convention)

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
