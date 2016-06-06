"""This module handles spectrum data objects."""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ..analysis.utils import resample

# STDLIB
import logging
logging.basicConfig(level=logging.INFO)
import re

# THIRD-PARTY
import numpy as np
from astropy.nddata import NDData, NDArithmeticMixin, NDIOMixin
from astropy.nddata.nduncertainty import StdDevUncertainty, NDUncertainty
from astropy.units import Unit, Quantity, spectral, spectral_density
from ..third_party.py_expression_eval import Parser


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

    # NOTE: Cannot have docstring here or Astropy will throw error!
    @classmethod
    def read(cls, *args, **kwargs):
        from ..interfaces.registries import io_registry

        return io_registry.read(cls, *args, **kwargs)

    def _from_self(self, other, copy_dispersion=False):
        """Create a new `Data` object using current property values."""
        return Data(name=self.name, data=other, unit=self.unit,
                    uncertainty=StdDevUncertainty(self.uncertainty),
                    mask=self.mask, wcs=self.wcs,
                    dispersion=self.dispersion if copy_dispersion else None,
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
    name : str
        Short description.
    """
    def __init__(self, source, mask=None, parent=None, name=''):
        self._source = source
        self._mask = mask.astype(bool) if mask is not None else np.ones(
            source.data.shape, dtype=bool)
        self._parent = parent
        self.name = self._source.name + " Layer" if not name else name
        self.units = (self._source.dispersion_unit,
                      self._source.unit if self._source.unit is not None
                      else Unit(""))

    @classmethod
    def new(cls, data, mask=None):
        new_data = Data(data)
        return Layer(new_data, mask)

    def _arithmetic(self, operator, other, propagate=True):
        if isinstance(other, Layer):
            # Make sure units are compatible
            if not other.data.unit.is_equivalent(
                    self.data.unit,
                    equivalencies=spectral_density(other.dispersion)):
                logging.error("Spectral data objects have incompatible units.")
                return

            # Create temporary arrays from the source object; this obviates
            # the need to re-create masks
            this_disp_arr = Quantity(self._source.dispersion,
                                     self._source.dispersion_unit)
            this_data_arr = Quantity(self._source.data, self._source.unit)

            if isinstance(self, ModelLayer):
                this_disp_arr[self._mask] = self.dispersion
                this_data_arr[self._mask] = self.data

            other_disp_arr = Quantity(other._source.dispersion,
                                      other._source.dispersion_unit).to(
                this_disp_arr.unit)
            other_data_arr = Quantity(other._source.data,
                                      other._source.unit).astype(
                other.data.dtype)

            if isinstance(other, ModelLayer):
                other_disp_arr[other._mask] = other.dispersion.to(
                    this_disp_arr.unit)
                other_data_arr[other._mask] = other.data

            disp_min = max(self.dispersion.value[0], other.dispersion.value[0])
            disp_max = min(self.dispersion.value[-1], other.dispersion.value[-1])

            this_mask = ((this_disp_arr.value >= disp_min) &
                         (this_disp_arr.value <= disp_max))
            other_mask = ((other_disp_arr.value >= disp_min) &
                          (other_disp_arr.value <= disp_max))

            # Check sampling; re-sample if incompatible
            this_step = self.dispersion[1] - self.dispersion[0]
            other_step = other.dispersion[1] - other.dispersion[0]

            if this_step != other_step:
                other_data_val = resample(other_data_arr.value[other_mask],
                                          other_disp_arr.value[other_mask],
                                          this_disp_arr.value[this_mask])
                temp_data_val = other_data_arr.value
                temp_data_val[other_mask] = other_data_val
                other_data_val = temp_data_val
                other_mask = this_mask
            else:
                other_data_val = other_data_arr.value

            this_data = self._source._from_self(this_data_arr.value)
            other_data = other._source._from_self(other_data_val,
                                                  copy_dispersion=True)
        # Assume that the operand is a single number
        else:
            if isinstance(other, Quantity):
                other = other.value

            this_data = self._source
            new_other = np.empty(shape=this_data.data.shape)
            new_other.fill(float(other))
            other_data = Data(new_other)
            other_data._dispersion_unit = this_data.dispersion_unit

            # If the operation is addition/subtraction, the units have to
            # match for arithmetic, otherwise unit should be nothing
            if operator in ['add', 'subtract']:
                other_data._unit = this_data.unit

            other_mask = self._mask

        # Perform arithmetic operation
        operator = getattr(this_data, operator)

        # Check for WCS incompatibility
        if this_data.wcs != other_data.wcs:
            logging.warning("WCS objects are not equivalent; overriding "
                            "wcs information on 'other'.")
            tmp_wcs = other_data._wcs
            other_data._wcs = this_data.wcs
            result_source = operator(other_data, propagate_uncertainties=propagate)
            other_data._wcs = tmp_wcs
        else:
            result_source = operator(other_data, propagate_uncertainties=propagate)

        # For cases where the dispersion cannot be recalculated, copy the
        # dispersion array
        result_source._dispersion = other_data.dispersion

        if result_source.dispersion_unit.is_unity():
            result_source._dispersion_unit = other_data.dispersion_unit

        # Create a layer from the source data object
        result_layer = Layer(result_source, other_mask, self._parent,
                             self.name)

        return result_layer

    @classmethod
    def from_formula(cls, formula, layers):
        if not formula:
            return

        layers = layers
        new_layer = cls._evaluate(layers, formula)

        if new_layer is None:
            return

        new_layer._window = None
        new_layer._parent = None
        new_layer.name = "Resultant"

        return new_layer

    @classmethod
    def _evaluate(cls, layers, formula):
        """
        Parse a string into an arithmetic expression.

        Parameters
        ----------
        layers : list
            List of `Layer` objects that correspond to the given variables.
        formula : str
            A string describing the arithmetic operations to perform.
        """
        parser = Parser()

        for layer in layers:
            formula = formula.replace(layer.name,
                                      layer.name.replace(" ", "_"))

        try:
            expr = parser.parse(formula)
        except Exception as e:
            logging.error(e)
            return

        # Extract variables
        vars = expr.variables()

        # List the models in the same order as the variables
        # sorted_layers = [next(l for v in vars for l in layers
        #                       if l.name.replace(" ", "_") == v)]
        # sorted_layers = [l for v in vars for l in layers
        #                  if l.name.replace(" ", "_") == v]
        sorted_layers = []

        for v in vars:
            for l in layers:
                if l.name.replace(" ", "_") == v:
                    sorted_layers.append(l)
                    break

        if len(sorted_layers) != len(vars):
            logging.error("Incorrect layer arithmetic formula: the number "
                          "of layers does not match the number of variables.")

        try:
            result = parser.evaluate(expr.simplify({}).toString(),
                                     dict(pair for pair in
                                          zip(vars, sorted_layers)))
        except Exception as e:
            logging.error("While evaluating formula: {}".format(e))
            return

        return result

    def __add__(self, other):
        new_layer = self._arithmetic("add", other)

        return new_layer

    def __sub__(self, other):
        new_layer = self._arithmetic("subtract", other)

        return new_layer

    def __mul__(self, other):
        new_layer = self._arithmetic("multiply", other, propagate=True)

        return new_layer

    def __truediv__(self, other):
        new_layer = self._arithmetic("divide", other, propagate=True)

        return new_layer

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

        # This is needed to catch improper mask usage by layers.
        if self._mask.shape != self._source.dispersion.shape:
            raise ValueError(
                'Mask shape mismatch, expect {0} but get {1}'.format(
                    self._source.dispersion.shape, self._mask.shape))

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

        # TODO: Is this too hacky?
        # Handle bad mask when fitting goes awry.
        if self._mask.shape != self._source.dispersion.shape:
            self._mask = np.ones(self._source.dispersion.shape, dtype=np.bool)

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

    @classmethod
    def from_formula(cls, formula, models):
        result_model = cls._evaluate(formula, models)
        return result_model

    @classmethod
    def _evaluate(cls, models, formula):
        try:
            parser = Parser()
            expr = parser.parse(formula)
        except:
            return

        # Extract variables
        vars = expr.variables()

        # List the models in the same order as the variables
        sorted_models = [m for v in vars for m in models if m.name == v]

        if len(sorted_models) > len(vars):
            logging.error("Incorrect model arithmetic formula: the number "
                          "of models does not match the number of variables.")
            return
        elif len(sorted_models) < len(vars):
            extras = [x for x in vars if x not in [y.name for y in
                                                  sorted_models]]

            for extra in extras:
                matches = re.findall('([\+\*\-\/]?\s?{})'.format(extra), formula)

                for match in matches:
                    formula = formula.replace(match, "")

            expr = parser.parse(formula)
            vars = expr.variables()

        result = parser.evaluate(expr.simplify({}).toString(),
                                 dict(pair for pair in
                                      zip(vars, sorted_models)))

        return result
