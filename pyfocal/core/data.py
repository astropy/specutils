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
    """
    def __init__(self, data, dispersion=None, name="", *args, **kwargs):
        super(Data, self).__init__(data, *args, **kwargs)
        self._dispersion = dispersion
        self.name = name or "New Data Object"
        self._layers = []

    @classmethod
    def read(cls, *args, **kwargs):
        from ..interfaces.registries import io_registry

        return io_registry.read(cls, *args, **kwargs)

    @property
    def dispersion(self):
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
        try:
            return self.wcs.wcs.cunit[0]
        except AttributeError:
            logging.warning("No dispersion unit information in WCS.")

        return Unit("")


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
    def __init__(self, source, mask, parent=None):
        super(Layer, self).__init__()
        self._source = source
        self._mask = mask
        self._parent = parent
        self._model = None
        self.name = self._source.name + " Layer"
        self.units = (self._source.dispersion_unit,
                      self._source.unit if self._source.unit is not None else "")

    @property
    def data(self):
        if self._model is not None:
            data = self._model(self.dispersion)
        else:
            data = self._source.data[self._mask]

        return Quantity(data, unit=self._source.unit).to(self.units[1])

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
        return self._source.uncertainty

    @property
    def mask(self):
        return self._source.mask[self._mask]

    @property
    def wcs(self):
        return self._source.wcs

    @property
    def meta(self):
        return self._source.meta

    def set_model(self, model):
        self._model = model
