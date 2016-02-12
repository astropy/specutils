from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from astropy.nddata import NDData, NDDataBase, NDArithmeticMixin, NDIOMixin
from .events import EventHook
import numpy as np
import logging
from astropy.units import Unit, Quantity


class Data(NDIOMixin, NDArithmeticMixin, NDData):
    """
    Class of the base data container for all data (of type
    :class:`numpy.ndarray`) that is passed around in Pyfocal. It inherits from
    :class:`astropy.nddata.NDData` and provides functionality for arithmetic
    operations, I/O, and slicing.

    Parameters
    ----------
    data : ndarray
        Flux values.

    dispersion : ndarray or `None`
        Dispersion values. If not given, this is calculated from WCS.

    dispersion_unit : `~astropy.units.Unit` or `None`
        Dispersion unit. If not given, this is obtained from WCS.

    name : str
        Short description of the spectrum.

    args : tuple
        Additional positional arguments.

    kwargs : dict
        Additional keyword arguments.

    Examples
    --------
    >>> d = Data.read(
    ...     'generic_spectra.fits', filter='Generic Fits (*.fits *.mits)')
    >>> d = Data.read(
    ...     'generic_spectra.txt', filter='ASCII (*.txt *.dat)')

    """
    def __init__(self, data, dispersion=None, dispersion_unit=None, name="",
                 *args, **kwargs):
        super(Data, self).__init__(data, *args, **kwargs)
        self._dispersion = dispersion
        self._dispersion_unit = dispersion_unit
        self.name = name or "New Data Object"
        self._layers = []

    # NOTE: Cannot have docstring here or Astropy will throw error!
    @classmethod
    def read(cls, *args, **kwargs):
        from ..interfaces.registries import io_registry

        return io_registry.read(cls, *args, **kwargs)

    @property
    def dispersion(self):
        """Dispersion values."""
        if self._dispersion is None:
            self._dispersion = np.arange(self.data.size)

            # try:
            crval = self.wcs.wcs.crval[0]
            cdelt = self.wcs.wcs.cdelt[0]
            end = self.data.shape[0] * cdelt + crval
            self._dispersion = np.arange(crval, end, cdelt)
            # except:
            #     logging.warning("Invalid FITS headers; constructing default "
            #                     "dispersion array.")

        return self._dispersion

    @property
    def dispersion_unit(self):
        """Unit of dispersion."""
        if self._dispersion_unit is None:
            try:
                self._dispersion_unit = self.wcs.wcs.cunit[0]
            except AttributeError:
                logging.warning("No dispersion unit information in WCS.")
                self._dispersion_unit = Unit("")

        return self._dispersion_unit


class Layer(object):
    """
    Base class to handle layers in Pyfocal.

    A layer is a "view" into a :class:`pyfocal.core.data.Data` object. It does
    not hold any data itself, but instead contains a special `mask` object
    and reference to the original data.

    Since :class:`pyfocal.core.data.Layer` inherits from
    :class:`astropy.nddata.NDDataBase` and provides the
    :class:`astropy.nddata.NDArithmeticMixin` mixin, it is also possible to
    do arithmetic operations on layers.
    """
    def __init__(self, source, mask, parent=None, window=None):
        super(Layer, self).__init__()
        self._source = source
        self._mask = mask
        self._parent = parent
        self._window = window
        self.name = self._source.name + " Layer"
        self.units = (self._source.dispersion_unit,
                      self._source.unit if self._source.unit is not None else "")

    @property
    def data(self):
        data = self._source.data[self._mask]

        return Quantity(data, unit=self._source.unit).to(self.units[1])

    @property
    def unit(self):
        return self._source.unit

    @property
    def dispersion(self):
        return Quantity(self._source.dispersion[self._mask],
                        unit=self._source.dispersion_unit).to(
                self.units[0])

    @dispersion.setter
    def dispersion(self, value, unit=""):
        self._source._dispersion = value

    @property
    def uncertainty(self):
        return self._source.uncertainty[self._mask]

    @property
    def mask(self):
        return self._source.mask[self._mask]

    @property
    def wcs(self):
        return self._source.wcs

    @property
    def meta(self):
        return self._source.meta


class ModelLayer(object):
    def __init__(self, source, model, parent=None):
        self._source = source
        self._model = model
        self._data = None
        self._window = self._source._window
        self.name = self._source.name + " Model"

    @property
    def data(self):
        if self._data is None:
            self._data = self._model(self._source.dispersion.value)

        return Quantity(self._data,
                        unit=self._source.unit).to(self._source.units[1])

    @property
    def dispersion(self):
        return self._source.dispersion

    @property
    def uncertainty(self):
        return None #self._source.uncertainty

    @property
    def mask(self):
        return self._source.mask

    @property
    def wcs(self):
        return self._source.wcs

    @property
    def meta(self):
        return self._source.meta

    @property
    def model(self):
        return self._model

    @model.setter
    def model(self, value):
        self._model = value
        self._data = self._model(self._source.dispersion.value)

    @property
    def layer(self):
        return self._source
