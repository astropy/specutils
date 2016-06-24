from specutils.spectrum1d import Spectrum1D

import numpy as np
import logging

from astropy.nddata import NDArithmeticMixin, NDSlicingMixin, NDIOMixin
import astropy.units as u


class GenericSpectrum1D(NDIOMixin, NDSlicingMixin, NDArithmeticMixin,
                        Spectrum1D):
    """
    A generalized spectrum data object that implements several `NDData`
    features such as IO, arithmetic, and slicing.
    """

    # The Spectrum1D wcs_attribute object is useless to use, let's avoid
    #  using it
    _wcs_attributes = {}

    def __init__(self, data, name="", dispersion=None, dispersion_unit=None,
                 *args, **kwargs):
        super(GenericSpectrum1D, self).__init__(flux=data, *args, **kwargs)
        self._dispersion = dispersion
        self._dispersion_unit = dispersion_unit
        self.name = name or "New Data Object"

    def from_self(self, **kwargs):
        """
        Create a new `GenericSpectrum1D` object using current property
        values. Takes all the arguments a Spectrum1D expects, arguments that
        are not included use this instance's values.
        """
        self_kwargs = {"data": self.data, "dispersion": self.dispersion,
                       "unit": self.unit, "wcs": self.wcs,
                       "uncertainty": self.uncertainty, "mask": self.mask,
                       "meta": self.meta}

        self_kwargs.update(kwargs)

        return self.__class__(**self_kwargs, copy=True)

    @classmethod
    def from_array(cls, data, *args, **kwargs):
        return cls(name="Array Data", data=data, **args, **kwargs)

    @property
    def dispersion(self):
        """Dispersion values."""
        if self._dispersion is None:
            self._dispersion = np.arange(self.data.size)

            try:
                crval = self.wcs.wcs.crval[0]

                # RuntimeWarning: cdelt will be ignored since cd is present
                try:
                    cdelt = self.wcs.wcs.cd[0][0]
                except:
                    cdelt = self.wcs.wcs.cdelt[0]

                end = self.data.shape[0] * cdelt + crval
                num = (end - crval) / cdelt

                # TODO: the values for the keywords are not guaranteed to be
                #  at the first index
                if hasattr(self.wcs.wcs, 'ctype') and "log" \
                        in self.wcs.wcs.ctype[-1].lower():
                    self._dispersion = np.logspace(crval, end, num)
                else:
                    self._dispersion = np.arange(crval, end, cdelt)
            except:
                logging.warning("Invalid FITS headers; constructing default "
                                "dispersion array.")

        return self._dispersion

    @property
    def dispersion_unit(self):
        """Unit of dispersion."""
        if self._dispersion_unit is None:
            try:
                self._dispersion_unit = self.wcs.wcs.cunit[0]
            except AttributeError:
                logging.warning("No dispersion unit information in WCS.")
                self._dispersion_unit = u.Unit("")

        return self._dispersion_unit


if __name__ == '__main__':
    pass