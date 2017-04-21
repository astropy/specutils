import logging

import numpy as np
from astropy.nddata import NDDataRef
from astropy.wcs import WCS, WCSSUB_SPECTRAL
from astropy.units import Unit, Quantity
import astropy.units.equivalencies as eq


__all__ = ['Spectrum1D']

class Spectrum1D(NDDataRef):
    def __init__(self, flux, dispersion=None, wcs=None, unit=None, disp_unit=None,
                 *args, **kwargs):
        if not isinstance(flux, Quantity):
            flux = Quantity(flux, unit=unit or "Jy")

        if dispersion is not None and not isinstance(dispersion, Quantity):
            dispersion = Quantity(dispersion, unit=disp_unit or "Angstrom")

        # If dispersion had not been defined, attempt to use the wcs
        # information, if it exists
        if dispersion is None and wcs is not None:
            if isinstance(wcs, WCS):
                # Try to reference the spectral axis
                wcs_spec = wcs.sub([WCSSUB_SPECTRAL])

                # Check to see if it actually is a real coordinate description
                if wcs_spec.naxis == 0:
                    # It's not real, so attempt to get the spectral axis by
                    # specifying axis by integer
                    wcs_spec = wcs.sub([wcs.naxis])

                # Construct the dispersion array
                dispersion = wcs_spec.all_pix2world(
                    np.arange(flux.shape[0]), 0)[0]

                # Try to get the dispersion unit information
                try:
                    disp_unit = wcs.wcs.cunit[0]
                except AttributeError:
                    logging.warning("No dispersion unit information in WCS.")
                    disp_unit = Unit("")

                dispersion = dispersion * disp_unit

            self._dispersion = dispersion

        super(Spectrum1D, self).__init__(data=flux.value, unit=flux.unit,
                                         wcs=wcs, *args, **kwargs)

    @property
    def dispersion(self):
        return self._dispersion

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

    def to_dispersion(self, unit, rest, type="doppler_relativistic"):
        """
        Converts the dispersion array to the given
        wavelength/velocity/frequency space unit given the rest value.

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
        self.dispersion : ~`astropy.units.Quantity`
            The converted dispersion array in the new dispersion space.
        """
        if Unit(unit).is_equivalent("cm/s") or \
                self.dispersion.unit.is_equivalent("cm/s"):
            equiv = getattr(eq, type)(rest)
        else:
            equiv = eq.spectral()

        new_data = self._dispersion.to(unit, equivalencies=equiv)

        self._dispersion = new_data

        return self.dispersion