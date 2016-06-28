from specutils.spectrum1d import Spectrum1D

import numpy as np
import logging

from astropy.nddata import NDArithmeticMixin, NDSlicingMixin, NDIOMixin
import astropy.units as u
from astropy.wcs import WCS


class GenericSpectrum1D(NDIOMixin, NDSlicingMixin, NDArithmeticMixin,
                        Spectrum1D):
    """
    A generalized spectrum data object that implements several `NDData`
    features such as IO, arithmetic, and slicing.
    """

    # The Spectrum1D wcs_attribute object is useless to us, let's avoid it
    _wcs_attributes = {}

    def __init__(self, data, name="", dispersion=None, dispersion_unit=None,
                 uncertainty=None, wcs=None, *args, **kwargs):
        super(GenericSpectrum1D, self).__init__(flux=data,
                                                wcs=wcs,
                                                uncertainty=uncertainty,
                                                indexer="UNUSED",
                                                *args, **kwargs)
        self._dispersion = dispersion
        self._dispersion_unit = dispersion_unit
        self.name = name or "New Data Object"

    def from_self(self, copy=True, **kwargs):
        """
        Create a new `GenericSpectrum1D` object using current property
        values. Takes all the arguments a Spectrum1D expects, arguments that
        are not included use this instance's values.
        """
        self_kwargs = {"data": self._data, "dispersion": self._dispersion,
                       "unit": self.unit, "wcs": self.wcs,
                       "uncertainty": self._uncertainty, "mask": self.mask,
                       "meta": self.meta}

        self_kwargs.update(kwargs)

        return self.__class__(**self_kwargs, copy=copy)

    @classmethod
    def from_array(cls, data, *args, **kwargs):
        return cls(name="Array Data", data=data, **args, **kwargs)

    @property
    def dispersion(self):
        """Dispersion values."""
        if self._dispersion is None:
            self._dispersion = np.arange(self.data.size)

            if isinstance(self.wcs, WCS):
                try:
                    crval = self.wcs.wcs.crval[0]

                    # RuntimeWarning: cdelt will be ignored since cd is present
                    try:
                        cdelt = self.wcs.wcs.cd[0][0]
                    except AttributeError:
                        cdelt = self.wcs.wcs.cdelt[0]

                    end = self.data.shape[0] * cdelt + crval
                    num = (end - crval) / cdelt

                    # TODO: the values for the keywords are not guaranteed to
                    # be at the first index
                    if hasattr(self.wcs.wcs, 'ctype') and "log" \
                            in self.wcs.wcs.ctype[-1].lower():
                        self._dispersion = np.logspace(crval, end, num)
                    else:
                        self._dispersion = np.arange(crval, end, cdelt)
                except AttributeError:
                    logging.warning("Invalid FITS headers; constructing "
                                    "default dispersion array.")

        return self._dispersion

    @property
    def dispersion_unit(self):
        """Unit of dispersion."""
        if self._dispersion_unit is None:
            try:
                self._dispersion_unit = self.wcs.wcs.cunit[0]
            except AttributeError:
                logging.warning("No dispersion unit information in WCS.")

                try:
                    self._dispersion_unit = u.Unit(
                        self.meta['header']['cunit'][0])
                except KeyError:
                    logging.warning("No dispersion unit information in meta.")

                    self._dispersion_unit = u.Unit("")

        return self._dispersion_unit

    @dispersion_unit.setter
    def dispersion_unit(self, value):
        self._dispersion_unit = value

    def __add__(self, other):
        return self.add(other, handle_meta='first_found')

    def __sub__(self, other):
        return self.subtract(other, handle_meta='first_found')

    def __mul__(self, other):
        return self.multiply(other, handle_meta='first_found')

    def __truediv__(self, other):
        return self.divide(other, handle_meta='first_found')

    def __len__(self):
        return len(self.data)
