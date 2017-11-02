from __future__ import division

import logging

import numpy as np
from astropy.nddata import NDDataRef
from astropy.wcs import WCS, WCSSUB_SPECTRAL
from astropy.units import Unit, Quantity, dimensionless_unscaled
from astropy import units as u
from astropy.utils.decorators import lazyproperty

from .spectrum_mixin import OneDSpectrumMixin

__all__ = ['Spectrum1D']


class Spectrum1D(OneDSpectrumMixin, NDDataRef):
    """
    Spectrum container for 1D spectral data.
    """
    def __init__(self, flux, spectral_axis=None, wcs=None, unit=None,
                 spectral_axis_unit=None, *args, **kwargs):

        if not isinstance(flux, Quantity):
            flux = Quantity(flux, unit=unit or "Jy")

        if spectral_axis is not None and not isinstance(spectral_axis, Quantity):
            spectral_axis = Quantity(spectral_axis, unit=spectral_axis_unit or u.AA)

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
                    if spectral_axis_unit == dimensionless_unscaled:
                        logging.warning("Spectral_axis unit is dimensionless. Pls set one.")
                except AttributeError:
                    logging.warning("No spectral_axis unit information in WCS.")
                    spectral_axis_unit = Unit("")

                spectral_axis = spectral_axis * spectral_axis_unit

                if wcs.wcs.restfrq != 0:
                    self._rest_value = wcs.wcs.restfrq * u.Hz
                elif wcs.wcs.restwav != 0:
                    self._rest_value = wcs.wcs.restwav * u.AA

        super(Spectrum1D, self).__init__(data=flux.value, unit=flux.unit,
                                         wcs=wcs, *args, **kwargs)



    @property
    def frequency(self):
        return self.spectral_axis.to(u.GHz, u.spectral())

    @property
    def wavelength(self):
        return self.spectral_axis.to(u.AA, u.spectral())

    @property
    def energy(self):
        return self.spectral_axis.to(u.eV, u.spectral())

    @lazyproperty
    def bin_edges(self):
        return self.wcs.bin_edges()

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


        equiv = getattr(eq, 'doppler_{0}'.format(
            self.velocity_convention))(self.rest_value)

        new_data = self.spectral_axis.to(u.km/u.s, equivalencies=equiv)

        return new_data


    def spectral_resolution(true_dispersion, delta_dispersion, axis=-1):
        """Evaluate the probability distribution of the spectral resolution.

        For example, to tabulate a binned resolution function at 6000A
        covering +/-10A in 0.2A steps::

        R = spectrum1d.spectral_resolution(
        ... 6000 * u.Angstrom, np.linspace(-10, 10, 51) * u.Angstrom)
        assert R.shape == (50,)
        assert np.allclose(R.sum(), 1.)

        To build a sparse resolution matrix for true wavelengths 4000-8000A
        in 0.1A steps::

        R = spectrum1d.spectral_resolution(
        ... np.linspace(4000, 8000, 40001)[:, np.newaxis] * u.Angstrom,
        ... np.linspace(-10, +10, 201) * u.Angstrom)
        assert R.shape == (40000, 200)
        assert np.allclose(R.sum(axis=1), 1.)

        Parameters
        ----------
        true_dispersion : ~`astropy.units.Quantity`
            True value(s) of dispersion for which the resolution should be
            evaluated.
        delta_dispersion : ~`astropy.units.Quantity`
            Array of (observed - true) dispersion bin edges to integrate the
            resolution probability density over.
        axis : int
            Which axis of ``delta_dispersion`` contains the strictly increasing
            dispersion values to interpret as bin edges.  The dimension of
            ``delta_dispersion`` along this axis must be at least two.

        Returns
        -------
        numpy array
            Array of dimensionless probabilities calculated as the integral of
            P(observed | true) over each bin in (observed - true). The output
            shape is the result of broadcasting the input shapes.
        """
        pass
