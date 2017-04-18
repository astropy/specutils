import logging

import numpy as np
from astropy.nddata import NDDataRef
from astropy.wcs import WCS, WCSSUB_SPECTRAL
from astropy.units import Unit


__all__ = ['Spectrum1D']

class Spectrum1D(NDDataRef):
    def __init__(self, flux, dispersion=None, *args, **kwargs):
        super(Spectrum1D, self).__init__(data=flux, *args, **kwargs)

    @property
    def flux(self):
        """
        Convenience property to access the underlying `data` array.

        Returns
        -------
        ndarray
            The spectrum data axis array.
        """
        return self.data

    @property
    def dispersion(self):
        """
        The dispersion axis of the spectral object. This property will return a
        default dispersion if no custom dispersion has been provided. The
        returned dispersion object will also be rebinned if a binning matrix
        has been calculated.

        Returns
        -------
        dispersion : ndarray
            The spectral dispersion values.
        """
        # If dispersion has not yet been defined, attempt to use the wcs
        # information, if it exists
        if self._dispersion is None and self.wcs is not None:
            self._dispersion = np.arange(self.data.shape[0])

            if isinstance(self.wcs, WCS):
                # Try to reference the spectral axis
                wcs_spec = self.wcs.sub([WCSSUB_SPECTRAL])

                # Check to see if it actually is a real coordinate description
                if wcs_spec.naxis == 0:
                    # It's not real, so attempt to get the spectral axis by
                    # specifying axis by integer
                    wcs_spec = self.wcs.sub([self.wcs.naxis])

                # Construct the dispersion array
                self._dispersion = wcs_spec.all_pix2world(
                    np.arange(self.data.shape[0]), 0)[0]

        return self._dispersion

    @property
    def dispersion_unit(self):
        # If wcs information is provided, attempt to get the dispersion unit
        # from the header
        if self._dispersion_unit is None and self.wcs is not None:
            try:
                self._dispersion_unit = self.wcs.wcs.cunit[0]
            except AttributeError:
                logging.warning("No dispersion unit information in WCS.")

                try:
                    self._dispersion_unit = Unit(
                        self.meta['header']['cunit'][0])
                except KeyError:
                    logging.warning("No dispersion unit information in meta.")

                    self._dispersion_unit = Unit("")

        return self._dispersion_unit