from copy import deepcopy
import warnings

import numpy as np
import astropy.units.equivalencies as eq
from astropy import units as u
from astropy.coordinates import SpectralCoord
from astropy.nddata import StdDevUncertainty
from astropy.utils.decorators import deprecated
from astropy.utils.exceptions import AstropyUserWarning
from astropy.wcs import WCS

from ..utils.wcs_utils import gwcs_from_array
from .spectral_axis import SpectralAxis, observer_for_frame
from .spectral_frame import SpectralMedium, SPECTRAL_FRAMES, normalize_frame

DOPPLER_CONVENTIONS = {}
DOPPLER_CONVENTIONS['radio'] = u.doppler_radio
DOPPLER_CONVENTIONS['optical'] = u.doppler_optical
DOPPLER_CONVENTIONS['relativistic'] = u.doppler_relativistic

__all__ = ['OneDSpectrumMixin', 'RedshiftMixin', 'SpectralFrameMixin']


class SpectralFrameMixin():
    '''
    Mixin exposing the medium, reference frame and observation metadata
    carried by the `~specutils.SpectralAxis` of a `~specutils.Spectrum` or
    `~specutils.SpectrumCollection`.
    '''

    @property
    def medium(self):
        """
        The `~specutils.spectra.spectral_frame.SpectralMedium` the wavelengths
        of the spectral axis are expressed in, or `None` if unknown.
        """
        return self.spectral_axis.medium

    @property
    def frame(self):
        """
        FITS ``SPECSYS`` code of the reference frame the spectral axis values
        are measured in (e.g. ``'TOPOCENT'``, ``'BARYCENT'``, ``'SOURCE'``),
        or `None` if unknown. See
        `~specutils.spectra.spectral_frame.SPECTRAL_FRAMES`.

        The ``radial_velocity`` and ``redshift`` of the spectrum are the
        velocity of the source relative to an observer at rest in this frame,
        i.e. the shift still to be applied to reach the ``'SOURCE'`` (rest)
        frame, where they are zero.
        """
        return self.spectral_axis.frame

    @property
    def in_rest_frame(self):
        """
        `True` if the spectral axis is in the rest frame of the source
        (``frame == 'SOURCE'``), `False` if it is in another frame, and `None`
        if the frame is unknown.
        """
        frame = self.frame
        return None if frame is None else frame == 'SOURCE'

    @property
    def target(self):
        """
        The position (and, if known, velocity) of the source as an
        `~astropy.coordinates.SkyCoord` or coordinate frame, or `None`.
        """
        return self.spectral_axis.target

    @property
    def observer(self):
        """
        The position and velocity of the observer whose rest frame the spectral
        axis is expressed in, or `None`.
        """
        return self.spectral_axis.observer

    @property
    def obstime(self):
        """The mid-point of the observation as a `~astropy.time.Time`, or `None`."""
        return self.spectral_axis.obstime

    @property
    def location(self):
        """
        Where the spectrum was recorded, as an
        `~astropy.coordinates.EarthLocation`, or `None`.
        """
        return self.spectral_axis.location

    @property
    def barycentric_correction(self):
        """
        The velocity to add to radial velocities measured in the current
        `frame` to obtain barycentric radial velocities, computed from the
        `observer` and `target`. This is the standard barycentric correction
        for a topocentric spectrum, and zero once the spectrum is in the
        ``'BARYCENT'`` frame. `None` if the observer or target is unknown.
        """
        axis = self.spectral_axis
        if axis.observer is None or axis.target is None:
            return None
        probe = SpectralCoord(1 * u.AA, observer=axis.observer, target=axis.target)
        shifted = probe.with_observer_stationary_relative_to('icrs')
        return (shifted.radial_velocity - probe.radial_velocity).to(u.km / u.s)

    def _with_spectral_axis(self, spectral_axis):
        """
        Return a copy of this object with ``spectral_axis`` in place of the
        current one. Implemented by the concrete classes.
        """
        raise NotImplementedError

    def with_medium(self, medium, scheme='inversion'):
        """
        Return a copy of this spectrum with the spectral axis converted to a
        different `~specutils.spectra.spectral_frame.SpectralMedium`.

        Parameters
        ----------
        medium : `~specutils.spectra.spectral_frame.SpectralMedium`, {'vacuum', 'air'} or dict
            The medium to convert to. For air, the refraction formula and air
            conditions of the target medium are used for the conversion; when
            converting from air, those of the current medium are.
        scheme : str, optional
            How to invert the refractive index when converting from air, see
            `~specutils.utils.wcs_utils.air_to_vac`.

        Returns
        -------
        `~specutils.Spectrum`
            A copy in the new medium. The original WCS is stored in
            ``meta['original_wcs']`` and replaced by a lookup table.
        """
        from ..utils.wcs_utils import vac_to_air, air_to_vac

        new_medium = SpectralMedium.from_input(medium)
        if new_medium is None:
            raise ValueError("A medium must be given.")
        current = self.medium
        if current is None:
            raise ValueError("The medium of this spectrum is unknown, so it cannot be "
                             "converted. Specify it with the ``medium`` argument when "
                             "creating the spectrum.")

        axis = self.spectral_axis
        if current == new_medium:
            return self._with_spectral_axis(axis)
        if axis.unit is u.pixel:
            raise u.UnitsError("Cannot convert the medium of a spectral axis in pixel units.")
        if not axis.unit.is_equivalent(u.m):
            raise u.UnitsError("The spectral axis must be in wavelength units to convert "
                               f"between media, not '{axis.unit}'. Use "
                               "``with_spectral_axis_unit`` first.")

        wavelengths = axis.quantity
        if current.is_air:
            wavelengths = air_to_vac(wavelengths, scheme=scheme, **current.refraction_kwargs)
        if new_medium.is_air:
            wavelengths = vac_to_air(wavelengths, **new_medium.refraction_kwargs)

        return self._with_spectral_axis(axis.replicate(value=wavelengths, medium=new_medium))

    def with_frame(self, frame, velocity=None):
        """
        Return a copy of this spectrum with the spectral axis transformed to a
        different reference frame.

        When both `observer` and `target` are known (for instance because the
        spectrum was created with ``frame``, ``location``, ``obstime`` and
        ``target``), the transformation is computed from them with
        `~astropy.coordinates.SpectralCoord.with_observer_stationary_relative_to`.
        Otherwise, a transformation to the ``'SOURCE'`` (rest) frame applies
        the spectrum's `radial_velocity`, and any other transformation
        requires ``velocity``.

        Parameters
        ----------
        frame : str
            FITS ``SPECSYS`` code of the frame to transform to, e.g.
            ``'BARYCENT'`` or ``'SOURCE'``. See
            `~specutils.spectra.spectral_frame.SPECTRAL_FRAMES`.
        velocity : `~astropy.units.Quantity` ['speed'], optional
            The correction to add to radial velocities measured in the current
            frame to obtain those in ``frame`` (for example the barycentric
            correction reported by a pipeline, when going from
            ``'TOPOCENT'`` to ``'BARYCENT'``). A positive value shifts the
            spectral axis to longer wavelengths.

        Returns
        -------
        `~specutils.Spectrum`
            A copy in the new frame. Its `radial_velocity` is the velocity of
            the source relative to the new frame (zero in ``'SOURCE'``). The
            original WCS is stored in ``meta['original_wcs']`` and replaced by
            a lookup table.
        """
        frame = normalize_frame(frame)
        if frame is None:
            raise ValueError("A frame must be given.")
        axis = self.spectral_axis
        if axis.unit is u.pixel:
            raise u.UnitsError("Cannot transform the frame of a spectral axis in pixel units.")
        if axis.medium is not None and axis.medium.is_air:
            warnings.warn("Applying a velocity shift to air wavelengths; the refractive "
                          "index is assumed constant over the shift. Convert to vacuum "
                          "first for full accuracy.", AstropyUserWarning)

        if frame == self.frame:
            if velocity is not None:
                raise ValueError(f"The spectrum is already in the '{frame}' frame; "
                                 "cannot also apply a velocity.")
            return self._with_spectral_axis(axis)

        has_observer_and_target = axis.observer is not None and axis.target is not None

        if velocity is not None:
            velocity = u.Quantity(velocity)
            if not velocity.unit.is_equivalent(u.km / u.s):
                raise u.UnitsError("velocity must have units of speed.")
            if has_observer_and_target:
                # Move the observer along the line of sight so that the
                # radial velocity changes by +velocity
                new_axis = axis.with_radial_velocity_shift(observer_shift=-velocity)
            else:
                new_axis = axis.with_radial_velocity_shift(target_shift=velocity)
                new_axis = self._rebuild_observer(new_axis, frame)
        elif has_observer_and_target:
            if frame == 'SOURCE':
                reference = axis.target
            elif frame == 'TOPOCENT':
                if axis.location is None or axis.obstime is None:
                    raise ValueError("Transforming to the 'TOPOCENT' frame requires the "
                                     "location and obstime of the observation.")
                reference = observer_for_frame(frame, axis.location, axis.obstime)
            else:
                reference = SPECTRAL_FRAMES[frame]
                if reference is None:
                    raise ValueError(f"astropy has no coordinate frame for '{frame}'; "
                                     "pass the ``velocity`` to apply instead.")
            new_axis = axis.with_observer_stationary_relative_to(reference)
        elif frame == 'SOURCE':
            if self.frame is None:
                raise ValueError("The current frame of the spectrum is unknown. Specify "
                                 "it with the ``frame`` argument when creating the spectrum.")
            new_axis = axis.to_rest()
        else:
            raise ValueError(
                f"Cannot transform from '{self.frame}' to '{frame}' without either the "
                "``velocity`` to apply, or both an observer and a target (create the "
                "spectrum with ``target``, ``location`` and ``obstime``).")

        return self._with_spectral_axis(new_axis.replicate(frame=frame))

    @staticmethod
    def _rebuild_observer(axis, frame):
        """
        After a manual velocity shift of an axis without a target, replace any
        observer with one at rest in ``frame`` so that the observer stays
        consistent with the frame, or drop it if that is not possible.
        """
        if axis.observer is None:
            return axis
        observer = None
        if axis.location is not None and axis.obstime is not None:
            observer = observer_for_frame(frame, axis.location, axis.obstime)
        if observer is None:
            warnings.warn(f"Dropping the observer of the spectral axis, which cannot be "
                          f"transformed to the '{frame}' frame.", AstropyUserWarning)
        return SpectralAxis(axis.quantity, radial_velocity=axis.radial_velocity,
                            doppler_rest=axis.doppler_rest,
                            doppler_convention=axis.doppler_convention,
                            observer=observer, **axis._metadata)

    def to_rest(self):
        """
        Return a copy of this spectrum in the rest frame of the source
        (``frame == 'SOURCE'``), with `radial_velocity` zero. Equivalent to
        ``with_frame('SOURCE')``. Unlike `shift_spectrum_to`, this records the
        frame of the result and does not modify the spectrum in place.
        """
        return self.with_frame('SOURCE')


class RedshiftMixin():
    '''
    Mixin to define properties and methods related to redshift and radial
    velocity that are common to `~specutils.Spectrum` and `~specutils.SpectrumCollection`.
    '''

    @property
    def spectral_axis(self):
        """
        Returns the SpectralCoord object.
        """
        return self._spectral_axis

    @property
    def redshift(self):
        """
        The redshift(s) of the objects represented by this spectrum.  May be
        scalar (if this spectrum's ``flux`` is 1D) or vector.  Note that
        the concept of "redshift of a spectrum" can be ambiguous, so the
        interpretation is set to some extent by either the user, or operations
        (like template fitting) that set this attribute when they are run on
        a spectrum.
        """
        return self.spectral_axis.redshift

    @property
    def radial_velocity(self):
        """
        The radial velocity(s) of the objects represented by this spectrum.  May
        be scalar (if this spectrum's ``flux`` is 1D) or vector.  Note that
        the concept of "RV of a spectrum" can be ambiguous, so the
        interpretation is set to some extent by either the user, or operations
        (like template fitting) that set this attribute when they are run on
        a spectrum.
        """
        return self.spectral_axis.radial_velocity

    def set_redshift_to(self, redshift):
        """
        This sets the redshift of the spectrum to be `redshift` *without*
        changing the values of the `spectral_axis`.

        If you want to shift the `spectral_axis` based on this value, use
        `shift_spectrum_to`.
        """
        new_spec_coord = self.spectral_axis.replicate(redshift=redshift)
        self._spectral_axis = new_spec_coord

    def set_radial_velocity_to(self, radial_velocity):
        """
        This sets the radial velocity of the spectrum to be `radial_velocity`
        *without* changing the values of the `spectral_axis`.

        If you want to shift the `spectral_axis` based on this value, use
        `shift_spectrum_to`.
        """
        new_spec_coord = self.spectral_axis.replicate(
            radial_velocity=radial_velocity
        )
        self._spectral_axis = new_spec_coord

    def shift_spectrum_to(self, *, redshift=None, radial_velocity=None):
        """
        This shifts in-place the values of the `spectral_axis`, given either a
        redshift or radial velocity.

        If you do *not* want to change the `spectral_axis`, use
        `set_redshift_to` or `set_radial_velocity_to`.
        """
        if redshift is not None and radial_velocity is not None:
            raise ValueError(
                "Only one of redshift or radial_velocity can be used."
            )

        old_redshift = self.redshift

        if redshift is not None:
            # with_radial_velocity_shift(redshift) looks wrong but astropy SpectralCoord handles
            # redshift input to that method
            new_spectral_axis = self.spectral_axis.with_radial_velocity_shift(
                -self.spectral_axis.radial_velocity
            ).with_radial_velocity_shift(redshift)
            self._spectral_axis = new_spectral_axis
        elif radial_velocity is not None:
            if not radial_velocity.unit.is_equivalent(u.km/u.s):
                raise u.UnitsError("Radial velocity must be a velocity.")

            new_spectral_axis = self.spectral_axis.with_radial_velocity_shift(
                -self.spectral_axis.radial_velocity
            ).with_radial_velocity_shift(radial_velocity)
            self._spectral_axis = new_spectral_axis
            redshift = radial_velocity.to(u.Unit(''), u.doppler_redshift())
        else:
            raise ValueError("One of redshift or radial_velocity must be set.")

        # Also store an updated WCS if we can update it.
        if isinstance(self.wcs, WCS):
            wcs_spectral_index = self.wcs.wcs.spec + 1
            h = self.wcs.to_header()
            spec_ctype = h[f'CTYPE{wcs_spectral_index}']
            z_factor = (1 + redshift) / (1 + old_redshift)
            if spec_ctype[0:4] != 'WAVE':
                # Frequency, wavenumber and energy all invert this factor. Note that the FITS
                # keyword for wavenumber is WAVN, which won't match here.
                z_factor = 1 / z_factor
            new_crval = h[f'CRVAL{wcs_spectral_index}'] * z_factor
            h[f'CRVAL{wcs_spectral_index}'] = new_crval.value
            pc_key = f'PC{wcs_spectral_index}_{wcs_spectral_index}'
            if pc_key in h:
                h[pc_key] *= z_factor
            if f'CDELT{wcs_spectral_index}' in h:
                new_cdelt = h[f'CDELT{wcs_spectral_index}'] * z_factor
                h[f'CDELT{wcs_spectral_index}'] = new_cdelt.value
            # WCS doesn't allow updating, but you can set it to None and then assign a new value
            self.wcs = None
            self.wcs = WCS(h, preserve_units=True)
        elif self.wcs is not None:
            # I don't know how to update a GWCS cleanly so for now, we replace it and store the
            # old one to retain any spatial information in the original
            self._original_wcs = self.wcs
            self.wcs = None
            self.wcs = gwcs_from_array(new_spectral_axis, self.flux.shape,
                                    spectral_axis_index=self.spectral_axis_index)


class OneDSpectrumMixin():
    @property
    def spectral_axis_index(self):
        return self._spectral_axis_index

    @property
    def _spectral_axis_len(self):
        """
        How many elements are in the spectral dimension?
        """
        return self.data.shape[self._spectral_axis_numpy_index]

    @property
    def _data_with_spectral_axis_last(self):
        """
        Returns a view of the data with the spectral axis last
        """
        if self._spectral_axis_numpy_index == self.data.ndim - 1:
            return self.data
        else:
            return self.data.swapaxes(self._spectral_axis_numpy_index,
                                      self.data.ndim - 1)

    @property
    def _data_with_spectral_axis_first(self):
        """
        Returns a view of the data with the spectral axis first
        """
        if self._spectral_axis_numpy_index == 0:
            return self.data
        else:
            return self.data.swapaxes(self._spectral_axis_numpy_index, 0)

    @property
    def spectral_wcs(self):
        """
        Returns the spectral axes of the WCS
        """
        return self.wcs.axes.spectral

    @property
    def spectral_axis(self):
        """
        Returns the SpectralCoord object.
        """
        return self._spectral_axis

    @property
    def rest_value(self):
        return self.spectral_axis.doppler_rest

    @rest_value.setter
    def rest_value(self, value):
        self.spectral_axis.doppler_rest = value

    @property
    def flux(self):
        """
        Converts the stored data and unit information into a quantity.

        Returns
        -------
        `~astropy.units.Quantity`
            Spectral data as a quantity.
        """
        return u.Quantity(self.data, unit=self.unit, copy=False)

    @deprecated('v1.13', alternative="with_flux_unit")
    def new_flux_unit(self, unit, equivalencies=None, suppress_conversion=False):
        return self.with_flux_unit(unit, equivalencies=equivalencies,
                                  suppress_conversion=suppress_conversion)

    def _convert_flux(self, unit, equivalencies=None, suppress_conversion=False):
        """This is always done in-place.
        Also see :meth:`with_flux_unit`."""

        if not suppress_conversion:
            if equivalencies is None:
                equivalencies = eq.spectral_density(self.spectral_axis)

            new_data = self.flux.to(unit, equivalencies=equivalencies)

            self._data = new_data.value
            self._unit = new_data.unit
        else:
            self._unit = u.Unit(unit)

        if self.uncertainty is not None:
            self.uncertainty = StdDevUncertainty(
                self.uncertainty.represent_as(StdDevUncertainty).quantity.to(
                    unit, equivalencies=equivalencies))

    def with_flux_unit(self, unit, equivalencies=None, suppress_conversion=False):
        """Returns a new spectrum with a different flux unit.
        If uncertainty is defined, it will be converted to
        `~astropy.nddata.StdDevUncertainty` in the new unit.

        Parameters
        ----------
        unit : str or `~astropy.units.Unit`
            The unit to convert the flux array to.

        equivalencies : list of equivalencies
            Custom equivalencies to apply to conversions.
            Set to spectral_density by default.

        suppress_conversion : bool
            Set to `True` if updating the flux unit without
            converting data values. This is ignored for
            ``uncertainty`` component.

        Returns
        -------
        new_spec : `~specutils.Spectrum`
            A new spectrum with the converted flux array
            (and uncertainty, if applicable).

        """
        new_spec = deepcopy(self)
        new_spec._convert_flux(
            unit, equivalencies=equivalencies, suppress_conversion=suppress_conversion)
        return new_spec

    @property
    def velocity_convention(self):
        """
        Returns the velocity convention
        """
        return self.spectral_axis.doppler_convention

    def with_velocity_convention(self, velocity_convention):
        new_spectral_axis = self.spectral_axis.replicate(
            doppler_convention=velocity_convention)
        return self.__class__(flux=self.flux, spectral_axis=new_spectral_axis, wcs=self.wcs,
                              meta=self.meta, uncertainty=self.uncertainty, mask=self.mask,
                              spectral_axis_index=self.spectral_axis_index)

    @property
    def rest_value(self):
        return self.spectral_axis.doppler_rest

    @rest_value.setter
    def rest_value(self, value):
        self.spectral_axis.doppler_rest = value

    @property
    def velocity(self):
        """
        Converts the spectral axis array to the given velocity space unit given
        the rest value.

        These aren't input parameters but required Spectrum attributes

        Parameters
        ----------
        unit : str or ~`astropy.units.Unit`
            The unit to convert the dispersion array to.
        rest : ~`astropy.units.Quantity`
            Any quantity supported by the standard spectral equivalencies
            (wavelength, energy, frequency, wave number).
        type : {"doppler_relativistic", "doppler_optical", "doppler_radio"}
            The type of doppler spectral equivalency.
        redshift or radial_velocity
            If present, this shift is applied to the final output velocity to
            get into the rest frame of the object.

        Returns
        -------
        new_data : `~astropy.units.Quantity`
            The converted dispersion array in the new dispersion space.
        """
        if self.rest_value is None:
            raise ValueError("Cannot get velocity representation of spectral "
                             "axis without specifying a reference value.")
        if self.velocity_convention is None:
            raise ValueError("Cannot get velocity representation of spectral "
                             "axis without specifying a velocity convention.")

        equiv = getattr(u.equivalencies, 'doppler_{0}'.format(
            self.velocity_convention))(self.rest_value)

        new_data = self.spectral_axis.to(u.km/u.s, equivalencies=equiv).quantity

        # if redshift/rv is present, apply it:
        if self.spectral_axis.radial_velocity is not None:
            new_data += self.spectral_axis.radial_velocity

        return new_data

    @deprecated('v1.13', alternative="with_spectral_axis_unit")
    def with_spectral_unit(self, unit, velocity_convention=None,
                           rest_value=None):
        self.with_spectral_axis_unit(unit, velocity_convention=velocity_convention,
                                     rest_value=rest_value)

    def with_spectral_axis_unit(self, unit, velocity_convention=None, rest_value=None):
        """
        Returns a new spectrum with a different spectral axis unit. Note that this creates a new
        object using the converted spectral axis and thus drops the original WCS, if it existed,
        replacing it with a lookup-table :class:`~gwcs.wcs.WCS` based on the new spectral axis. The
        original WCS will be stored in the ``original_wcs`` entry of the new object's ``meta``
        dictionary.

        Parameters
        ----------
        unit : :class:`~astropy.units.Unit`
            Any valid spectral unit: velocity, (wave)length, or frequency.
            Only vacuum units are supported.
        velocity_convention : 'relativistic', 'radio', or 'optical'
            The velocity convention to use for the output velocity axis.
            Required if the output type is velocity. This can be either one
            of the above strings, or an `astropy.units` equivalency.
        rest_value : :class:`~astropy.units.Quantity`
            A rest wavelength or frequency with appropriate units.  Required if
            output type is velocity.  The spectrum's WCS should include this
            already if the *input* type is velocity, but the WCS's rest
            wavelength/frequency can be overridden with this parameter.

            .. note: This must be the rest frequency/wavelength *in vacuum*,
                     even if your spectrum has air wavelength units

        """
        velocity_convention = velocity_convention if velocity_convention is not None else self.velocity_convention  # noqa
        rest_value = rest_value if rest_value is not None else self.rest_value
        unit = self._new_wcs_argument_validation(unit, velocity_convention, rest_value)

        # Store the original unit information and WCS for posterity
        meta = deepcopy(self._meta)

        if 'original_spectral_axis_unit' not in self._meta:
            orig_unit = self.wcs.unit[0] if hasattr(self.wcs, 'unit') else self.spectral_axis.unit
            meta['original_spectral_axis_unit'] = orig_unit

        if 'original_wcs' not in self.meta:
            meta['original_wcs'] = self.wcs.deepcopy()

        new_spectral_axis = self.spectral_axis.to(unit, doppler_convention=velocity_convention,
                                                  doppler_rest=rest_value)

        return self.__class__(flux=self.flux, spectral_axis=new_spectral_axis, meta=meta,
                              uncertainty=self.uncertainty, mask=self.mask)

    def with_spectral_axis_and_flux_units(self, spectral_axis_unit, flux_unit,
                                          velocity_convention=None, rest_value=None,
                                          flux_equivalencies=None, suppress_flux_conversion=False):
        """Perform :meth:`with_spectral_axis_unit` and :meth:`with_flux_unit` together.
        See the respective methods for input and output definitions.

        Returns
        -------
        new_spec : `~specutils.Spectrum`
            Spectrum in requested units.

        """
        new_spec = self.with_spectral_axis_unit(
            spectral_axis_unit, velocity_convention=velocity_convention, rest_value=rest_value)
        new_spec._convert_flux(
            flux_unit, equivalencies=flux_equivalencies, suppress_conversion=suppress_flux_conversion)
        return new_spec

    def _axis_length_validation(self):
        pass

    def _new_wcs_argument_validation(self, unit, velocity_convention,
                                     rest_value):
        # Allow string specification of units, for example
        if not isinstance(unit, u.UnitBase):
            unit = u.Unit(unit)

        # Velocity conventions: required for frq <-> velo
        # convert_spectral_axis will handle the case of no velocity
        # convention specified & one is required
        if velocity_convention in DOPPLER_CONVENTIONS:
            velocity_convention = DOPPLER_CONVENTIONS[velocity_convention]
        elif (velocity_convention is not None and
              velocity_convention not in DOPPLER_CONVENTIONS.values()):
            raise ValueError("Velocity convention must be radio, optical, "
                             "or relativistic.")

        # If rest value is specified, it must be a quantity
        if (rest_value is not None and
            (not hasattr(rest_value, 'unit') or
             not rest_value.unit.is_equivalent(u.m, u.spectral()))):
            raise ValueError("Rest value must be specified as an astropy "
                             "quantity with spectral equivalence.")

        return unit

    def _check_strictly_increasing_decreasing(self):
        """
        Check that the self._spectral_axis is strictly increasing or decreasing
        and raise an error if its not.

        """

        spec_axis = self._spectral_axis

        sorted_increasing = np.all(spec_axis[1:] >= spec_axis[:-1])
        if sorted_increasing:  # check increasing first, probably most common case
            self._spectral_axis_direction = 'increasing'
            return True
        sorted_decreasing = np.all(spec_axis[1:] <= spec_axis[:-1])
        if sorted_decreasing:
            self._spectral_axis_direction = 'decreasing'
            return True
        return False


class InplaceModificationMixin:
    # Example methods follow to demonstrate how methods can be written to be
    # agnostic of the non-spectral dimensions.

    def substract_background(self, background):
        """
        Proof of concept, this subtracts a background spectrum-wise
        """

        data = self._data_with_spectral_axis_last

        if callable(background):
            # create substractable array
            pass
        elif (isinstance(background, np.ndarray) and
              background.shape == data[-1].shape):
            substractable_continuum = background
        else:
            raise ValueError(
                "background needs to be callable or have the same shape as the spectum")

        data[-1] -= substractable_continuum

    def normalize(self):
        """
        Proof of concept, this normalizes each spectral dimension based
        on a trapezoidal integration.
        """

        # this gets a view - if we want normalize to return a new NDData object
        # then we should make _data_with_spectral_axis_first return a copy.
        data = self._data_with_spectral_axis_first

        dx = np.diff(self.spectral_axis)
        dy = 0.5 * (data[:-1] + data[1:])

        norm = np.sum(dx * dy.transpose(), axis=-1).transpose()

        data /= norm

    def spectral_interpolation(self, spectral_value, flux_unit=None):
        """
        Proof of concept, this interpolates along the spectral dimension
        """

        data = self._data_with_spectral_axis_last

        from scipy.interpolate import interp1d

        interp = interp1d(self.spectral_axis.value, data)

        x = spectral_value.to(self.spectral_axis.unit,
                              equivalencies=u.spectral())
        y = interp(x)

        if self.unit is not None:
            y *= self.unit

        if flux_unit is None:  # Lim: Is this acceptable?
            return y
        else:
            return y.to(flux_unit, equivalencies=u.spectral_density(x))
