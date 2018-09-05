import logging
from copy import deepcopy

import numpy as np
from astropy.wcs import WCSSUB_SPECTRAL
from astropy.units import Unit
from astropy import units as u
import astropy.units.equivalencies as eq
from astropy.utils.decorators import lazyproperty

DOPPLER_CONVENTIONS = {}
DOPPLER_CONVENTIONS['radio'] = u.doppler_radio
DOPPLER_CONVENTIONS['optical'] = u.doppler_optical
DOPPLER_CONVENTIONS['relativistic'] = u.doppler_relativistic

__all__ = ['OneDSpectrumMixin']


class OneDSpectrumMixin(object):
    @property
    def _spectral_axis_numpy_index(self):
        return self.data.ndim - 1 - self.wcs.wcs.spec

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
        return self.wcs.axes.spectral

    @lazyproperty
    def spectral_axis(self):
        """
        Returns a Quantity array with the values of the spectral axis.
        """
        return self.wcs.pixel_to_world(np.arange(self.flux.shape[-1]))

    @property
    def spectral_axis_unit(self):
        """
        Returns the units of the spectral axis.
        """
        return self.wcs.spectral_axis_unit

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

    def new_flux_unit(self, unit, equivalencies=None, suppress_conversion=False):
        """
        Converts the flux data to the specified unit.  This is an in-place
        change to the object.

        Parameters
        ----------
        unit : str or `~astropy.units.Unit`
            The unit to conver the flux array to.

        equivalencies : list of equivalencies
            Custom equivalencies to apply to conversions.
            Set to spectral_density by default.

        suppress_conversion : bool
            Set to true if updating the unit without
            converting data values.

        Returns
        -------
        `~specutils.Spectrum1D`
            A new spectrum with the converted flux array
        """
        new_spec = deepcopy(self)
        if not suppress_conversion:
            if equivalencies is None:
                equivalencies = eq.spectral_density(self.spectral_axis)

            new_data = self.flux.to(
                unit, equivalencies=equivalencies)

            new_spec._data = new_data.value
            new_spec._unit = new_data.unit
        else:
            new_spec._unit = u.Unit(unit)

        return new_spec

    @property
    def velocity_convention(self):
        return self._velocity_convention

    def with_velocity_convention(self, velocity_convention):
        return self.__class__(flux=self.flux, wcs=self.wcs, meta=self.meta,
                              velocity_convention=velocity_convention)

    @property
    def rest_value(self):
        return self._rest_value

    @rest_value.setter
    def rest_value(self, value):
        if not hasattr(value, 'unit') or not value.unit.is_equivalent(u.Hz, u.spectral()):
            raise u.UnitsError(
                "Rest value must be energy/wavelength/frequency equivalent.")

        self._rest_value = value

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

        Returns
        -------
        ~`astropy.units.Quantity`
            The converted dispersion array in the new dispersion space.
        """
        if not hasattr(self, '_rest_value') or self._rest_value is None:
            raise ValueError("Cannot get velocity representation of spectral "
                             "axis without specifying a reference value.")
        if not hasattr(self, '_velocity_convention') or self._velocity_convention is None:
            raise ValueError("Cannot get velocity representation of spectral "
                             "axis without specifying a velocity convention.")

        equiv = getattr(u.equivalencies, 'doppler_{0}'.format(
            self.velocity_convention))(self.rest_value)

        new_data = self.spectral_axis.to(u.km/u.s, equivalencies=equiv)

        return new_data

    def with_spectral_unit(self, unit, velocity_convention=None,
                           rest_value=None):
        """
        Returns a new spectrum with a different spectral axis unit.

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
        new_wcs, new_meta = self._new_spectral_wcs(
            unit=unit,
            velocity_convention=velocity_convention or self._velocity_convention,
            rest_value=rest_value or self.rest_value)

        spectrum = self.__class__(flux=self.flux, wcs=new_wcs, meta=new_meta)

        return spectrum

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

    def _new_spectral_wcs(self, unit, velocity_convention=None,
                          rest_value=None):
        """
        Returns a new WCS with a different Spectral Axis unit.

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
            output type is velocity.  The cube's WCS should include this
            already if the *input* type is velocity, but the WCS's rest
            wavelength/frequency can be overridden with this parameter.

            .. note: This must be the rest frequency/wavelength *in vacuum*,
                     even if your cube has air wavelength units
        """

        unit = self._new_wcs_argument_validation(unit, velocity_convention,
                                                 rest_value)

        if velocity_convention is not None:
            equiv = getattr(u, 'doppler_{0}'.format(velocity_convention))
            rest_value.to(unit, equivalencies=equiv)

        # Store the original unit information for posterity
        meta = self._meta.copy()

        if 'original_unit' not in self._meta:
            meta['original_unit'] = self.wcs.spectral_axis_unit

        # Create the new wcs object
        new_wcs = self.wcs.with_spectral_unit(unit=unit,
                                         rest_value=rest_value,
                                         velocity_convention=velocity_convention)

        return new_wcs, meta


class InplaceModificationMixin(object):
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
