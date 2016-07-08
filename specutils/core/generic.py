from specutils.spectrum1d import Spectrum1D

import numpy as np
import logging

from astropy.nddata import NDArithmeticMixin, NDSlicingMixin, NDIOMixin
import astropy.units as u
from astropy.wcs import WCS, WCSSUB_SPECTRAL


class Spectrum1DRef(NDIOMixin, NDSlicingMixin, NDArithmeticMixin, Spectrum1D):
    """
    A generalized spectrum data object that implements several `NDData`
    features such as IO, arithmetic, and slicing.

    - IO is incorporated by the addition of `read` and `write` class
      methods that utilize the `~astropy.io.registry`. By default,
      it recognizes built-in astropy formats listed in the `unified read/write
      interface documentation <http://docs.astropy.org/en/stable/io/unified
      .html#table-io>`_.

    - The Arithmetic mixin on :class:`~astropy.nddata.NDData` objects allows
      for operations on the base `data` object as well as the `uncertainty`
      value, including error propagation when using an
      :class:`~astropy.nddata.NDUncertainty`-derived subclass. The
      mixin also handles the `WCS`, `meta` data, and other aspects of the
      :class:`~astropy.nddata.NDData` object.

    - Direct slicing of the :class:`~astropy.nddata.NDData` object is possible
      via the slicing mixin.

    This spectrum object differs from the :class:`~specutils.Spectrum1D` by
    relying on the :class:`~astropy.wcs.WCS` object to maintain WCS
    information in addition to the above listed extensions.
    """

    # The Spectrum1D wcs_attribute object is useless to us, let's avoid it
    _wcs_attributes = {}

    def __init__(self, data, name="", dispersion=None, dispersion_unit=None,
                 uncertainty=None, wcs=None, *args, **kwargs):
        super(Spectrum1DRef, self).__init__(flux=data,
                                            wcs=wcs,
                                            uncertainty=uncertainty,
                                            indexer="UNUSED",
                                            *args, **kwargs)
        self._dispersion = dispersion
        self._dispersion_unit = dispersion_unit
        self.name = name or "New Data Object"

    @classmethod
    def copy(cls, original, deep_copy=True, **kwargs):
        """
        Create a new `GenericSpectrum1D` object using current property
        values. Takes all the arguments a Spectrum1D expects, arguments that
        are not included use this instance's values.
        """
        self_kwargs = {"data": original._data,
                       "dispersion": original._dispersion,
                       "unit": original.unit, "wcs": original.wcs,
                       "uncertainty": original._uncertainty,
                       "mask": original.mask, "meta": original.meta}

        self_kwargs.update(kwargs)

        return cls(copy=deep_copy, **self_kwargs)

    @classmethod
    def from_array(cls, data, *args, **kwargs):
        return cls(name="Array Data", data=data, *args, **kwargs)

    @property
    def dispersion(self):
        """
        Override parent method in order to construct dispersion array from
        Astropy's WCS object instead of the custom specutils WCS object.
        """
        if self._dispersion is None:
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
