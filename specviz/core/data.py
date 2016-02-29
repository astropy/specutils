"""This module handles spectrum data objects."""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# STDLIB
import logging
logging.basicConfig(level=logging.INFO)
import numbers

# THIRD-PARTY
import numpy as np
from astropy.nddata import NDData, NDArithmeticMixin, NDIOMixin
from astropy.nddata.nduncertainty import StdDevUncertainty, NDUncertainty
from astropy.units import Unit, Quantity


class Data(NDIOMixin, NDArithmeticMixin, NDData):
    """Class of the base data container for all data (of type
    :class:`numpy.ndarray`) that is passed around in SpecViz. It inherits from
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
        super(Data, self).__init__(data=data, *args, **kwargs)

        # TODO: This is a temporary workaround until astropy releases its
        # next version which introduces changes to the NDData object
        if self.uncertainty is not None:
            self.uncertainty.parent_nddata = self

        self._dispersion = dispersion
        self._dispersion_unit = dispersion_unit
        self.name = name or "New Data Object"
        self._layers = []

    # NOTE: Cannot have docstring here or Astropy will throw error!
    @classmethod
    def read(cls, *args, **kwargs):
        from ..interfaces.registries import io_registry

        return io_registry.read(cls, *args, **kwargs)

    def _from_self(self, other):
        """Create a new `Data` object using current property values."""
        return Data(name=self.name, data=other, unit=self.unit,
                    uncertainty=StdDevUncertainty(self.uncertainty),
                    mask=self.mask, wcs=self.wcs, dispersion=self.dispersion,
                    dispersion_unit=self.dispersion_unit)

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
                self._dispersion_unit = Unit("")

        return self._dispersion_unit


class Layer(object):
    """Class to handle layers in SpecViz.

    A layer is a "view" into a :class:`Data` object. It does
    not hold any data itself, but instead contains a special ``mask`` object
    and reference to the original data.

    Since :class:`Data` inherits from
    :class:`astropy.nddata.NDDataBase` and provides the
    :class:`astropy.nddata.NDArithmeticMixin` mixin, it is also possible to
    do arithmetic operations on layers.

    Parameters
    ----------
    source : `Data`
        Spectrum data object.
    mask : ndarray
        Mask for the spectrum data.
    parent : obj or `None`
        GUI parent.
    window : obj or `None`
        GUI window.
    name : str
        Short description.
    """
    def __init__(self, source, mask, parent=None, name=''):
        self._source = source
        self._mask = mask
        self._parent = parent
        self.name = self._source.name + " Layer" if not name else name
        self.units = (self._source.dispersion_unit,
                      self._source.unit if self._source.unit is not None
                      else Unit(""))

    def _arithmetic(self, operator, other, propagate=True):
        # The operand is a single number
        if isinstance(other, numbers.Number):
            new = np.empty(shape=self.data.shape)
            new.fill(other)
            other = new
        # The operand is an array
        elif isinstance(other, np.ndarray) or isinstance(other, list):
            if isinstance(other, Quantity):
                other = other.value

            other = self._source._from_self(other)
        elif isinstance(other, Layer) or isinstance(other, ModelLayer):
            other = other._source

        if isinstance(other, Data):
            if self._source.wcs != other.wcs:
                logging.warning("WCS objects are not equivalent; overriding "
                                "wcs information on 'other'.".format())
                tmp_wcs = other._wcs
                other._wcs = self._source.wcs
                new_source = operator(other, propagate_uncertainties=propagate)
                other._wcs = tmp_wcs
            else:
                new_source = operator(other, propagate_uncertainties=propagate)

            # Force the same dispersion
            new_source._dispersion = np.copy(self._source._dispersion)

            # Apply arithmetic to the dispersion unit
            if operator.__name__ == 'multiply':
                new_source._dispersion_unit = self._source.dispersion_unit * \
                                              other.dispersion_unit
            elif operator.__name__ == 'divide':
                new_source._dispersion_unit = self._source.dispersion_unit / \
                                              other.dispersion_unit
            else:
                new_source._dispersion_unit = self._source._dispersion_unit
        else:
            new_source = operator(other, propagate_uncertainties=propagate)

        return new_source

    def __add__(self, other):
        new_source = self._arithmetic(self._source.add, other)

        return Layer(new_source, self._mask, self._parent, self.name)

    def __sub__(self, other):
        new_source = self._arithmetic(self._source.subtract, other)

        return Layer(new_source, self._mask, self._parent, self.name)

    def __mul__(self, other):
        new_source = self._arithmetic(self._source.multiply, other,
                                      propagate=True)

        return Layer(new_source, self._mask, self._parent, self.name)

    def __truediv__(self, other):
        new_source = self._arithmetic(self._source.divide, other,
                                      propagate=False)

        return Layer(new_source, self._mask, self._parent, self.name)

    @property
    def data(self):
        """Flux quantity with mask applied."""
        return Quantity(self._source.data[self._mask],
                        unit=self._source.unit).to(self.units[1])

    @property
    def unit(self):
        """Flux unit."""
        return self._source.unit

    @property
    def dispersion(self):
        """Dispersion quantity with mask applied."""
        return Quantity(self._source.dispersion[self._mask],
                        unit=self._source.dispersion_unit).to(self.units[0])

    @dispersion.setter
    def dispersion(self, value, unit=""):
        self._source._dispersion = value

    @property
    def uncertainty(self):
        """Flux uncertainty with mask applied."""
        return self._source.uncertainty[self._mask]

    @property
    def mask(self):
        """Mask for spectrum data."""
        return self._source.mask[self._mask]

    @property
    def wcs(self):
        """WCS for spectrum data."""
        return self._source.wcs

    @property
    def meta(self):
        """Spectrum metadata."""
        return self._source.meta


class ModelLayer(Layer):
    """A layer for spectrum with a model applied.

    Parameters
    ----------
    model : obj
        Astropy model.
    source : `Data`
        Spectrum data object.
    mask : ndarray
        Mask for the spectrum data.
    parent : obj or `None`
        GUI parent.
    window : obj or `None`
        GUI window.
    name : str
        Short description.
    """
    def __init__(self, model, source, mask, parent=None, name=''):
        name = source.name + " Model Layer" if not name else name
        super(ModelLayer, self).__init__(source, mask, parent, name)

        self._data = None
        self._model = model

        logging.info('Created ModelLayer object: {0}'.format(name))

    @property
    def data(self):
        """Flux quantity from model."""
        self._data = self._model(self.dispersion.value)

        return Quantity(self._data,
                        unit=self._source.unit).to(self.units[1])

    @property
    def uncertainty(self):
        """
        Models do not need to contain uncertainties; override parent
        class method.
        """
        return None

    @property
    def model(self):
        """Spectrum model."""
        return self._model

    @model.setter
    def model(self, value):
        self._model = value
        self._data = self._model(self.dispersion.value)
