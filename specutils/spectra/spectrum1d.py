import logging

import numpy as np
from astropy.nddata import NDDataRef
from astropy.wcs import WCS, WCSSUB_SPECTRAL
from astropy.units import Unit, Quantity
from astropy import units as u
import astropy.units.equivalencies as eq


__all__ = ['Spectrum1D']

class Spectrum1D(NDDataRef):
    """
    Spectrum container for 1D spectral data.
    """
    def __init__(self, flux, spectral_axis=None, wcs=None, unit=None, spectral_axis_unit=None,
                 *args, **kwargs):
        if not isinstance(flux, Quantity):
            flux = Quantity(flux, unit=unit or "Jy")

        if spectral_axis is not None and not isinstance(spectral_axis, Quantity):
            spectral_axis = Quantity(spectral_axis, unit=spectral_axis_unit or "Angstrom")

        # If spectral_axis had not been defined, attempt to use the wcs
        # information, if it exists
        if spectral_axis is None and wcs is not None:
            if isinstance(wcs, WCS):
                # Try to reference the spectral axis
                wcs_spec = wcs.sub([WCSSUB_SPECTRAL])

                # Check to see if it actually is a real coordinate description
                if wcs_spec.naxis == 0:
                    # It's not real, so attempt to get the spectral axis by
                    # specifying axis by integer
                    wcs_spec = wcs.sub([wcs.naxis])

                # Construct the spectral_axis array
                spectral_axis = wcs_spec.all_pix2world(
                    np.arange(flux.shape[0]), 0)[0]

                # Try to get the spectral_axis unit information
                try:
                    spectral_axis_unit = wcs.wcs.cunit[0]
                except AttributeError:
                    logging.warning("No spectral_axis unit information in WCS.")
                    spectral_axis_unit = Unit("")

                spectral_axis = spectral_axis * spectral_axis_unit

                if wcs.wcs.restfrq != 0:
                    self.rest_value = wcs.wcs.restfrq * u.Hz
                if wcs.wcs.restwav != 0:
                    self.rest_value = wcs.wcs.restwav * u.AA

        self._spectral_axis = spectral_axis

        super(Spectrum1D, self).__init__(data=flux.value, unit=flux.unit,
                                         wcs=wcs, *args, **kwargs)

    @property
    def spectral_axis(self):
        """
        The spectral axis data.

        Returns
        -------
        ~`astropy.units.Quantity`
            Spectral axis data as a quantity.
        """
        return self._spectral_axis

    @property
    def flux(self):
        """
        Converts the stored data and unit information into a quantity.

        Returns
        -------
        ~`astropy.units.Quantity`
            Spectral data as a quantity.
        """
        return self.data * Unit(self.unit)

    def to_flux(self, unit):
        """
        Converts the flux data to the specified unit.

        Parameters
        ----------
        unit : str or ~`astropy.units.Unit`
            The unit to conver the flux array to.

        Returns
        -------
        ~`astropy.units.Quantity`
            The converted flux array.
        """
        new_data = self.flux.to(
            unit, equivalencies=eq.spectral_density(self.spectral_axis))

        self._data = new_data.value
        self._unit = new_data.unit

        return self.flux

    @property
    def frequency(self):
        return self._spectral_axis.to(u.GHz, u.spectral())

    @property
    def wavelength(self):
        return self._spectral_axis.to(u.AA, u.spectral())

    @property
    def energy(self):
        return self._spectral_axis.to(u.ev, u.spectral())

    @property
    def velocity_convention(self):
        return self._velocity_convention

    @velocity_convention.setter
    def velocity_convention(self, value):
        if value not in ('relativistic', 'optical', 'radio'):
            raise ValueError("The allowed velocity conveintions are 'optical' "
                             "(linear with respect to wavelength), 'radio' "
                             "(linear with respect to frequency), and 'relativistic'.")
        self._velocity_convention = value

    @property
    def rest_value(self):
        return self._rest_value

    @rest_value.setter
    def rest_value(self, value):
        if not hasattr(value, 'unit') or not value.unit.is_equivalent(u.Hz, u.spectral()):
            raise ValueError("Rest value must be energy/wavelength/frequency equivalent.")
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
        if not hasattr(self, '_rest_value'):
            raise ValueError("Cannot get velocity representation of spectral "
                             "axis without specifying a reference value.")
        if not hasattr(self, '_velocity_convention'):
            raise ValueError("Cannot get velocity representation of spectral "
                             "axis without specifying a velocity convention.")


        equiv = getattr(eq, self.velocity_convention)('doppler_{0}'.format(self.rest_value))

        new_data = self.spectral_axis.to(u.km/u.s, equivalencies=equiv)

        return new_data
