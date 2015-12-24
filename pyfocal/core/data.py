from __future__ import absolute_import, division, print_function

from astropy.nddata import NDData, NDDataBase, NDArithmeticMixin, NDIOMixin
from .events import EventHook


class Data(NDIOMixin, NDArithmeticMixin, NDData):
    """
    Class of the base data container for all data (of type
    :class:`numpy.ndarray`) that is passed around in Pyfocal. It inherits from
    :class:`astropy.nddata.NDData` and provides functionality for arithmetic
    operations, I/O, and slicing.
    """
    def __init__(self, *args, **kwargs):
        super(Data, self).__init__(*args, **kwargs)
        self.name = "New Data Object"
        self._layers = []

    @classmethod
    def read(cls, *args, **kwargs):
        from ..interfaces.registries import io_registry

        return io_registry.read(cls, *args, **kwargs)


class Layer(NDArithmeticMixin, NDDataBase):
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
        self.name = self._source.name + " Layer"

    @property
    def data(self):
        return self._source.data[self._mask]

    @property
    def mask(self):
        return self._source.mask[self._mask]

    @property
    def unit(self):
        return self._source.unit

    @property
    def wcs(self):
        return self._source.wcs

    @property
    def meta(self):
        return self._source.meta
