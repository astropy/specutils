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

        self._dispersion = dispersion
        self._dispersion_unit = dispersion_unit
        self.name = name or "New Data Object"

    # NOTE: Cannot have docstring here or Astropy will throw error!
    @classmethod
    def read(cls, *args, **kwargs):
        from ..interfaces.registries import io_registry

        return io_registry.read(cls, *args, **kwargs)

    def copy(self, data, dispersion=None, dispersion_unit=None, mask=None,
             unit=None, uncertainty=None):
        """Create a new `Data` object using current property values."""
        return Data(name=self.name, data=data,
                    unit=unit if unit is not None else self.unit,
                    uncertainty=StdDevUncertainty(
                        uncertainty if uncertainty is not None else
                        self.uncertainty),
                    mask=mask if mask is not None else self.mask,
                    wcs=self.wcs,
                    dispersion=dispersion if dispersion is not None else
                                             self.dispersion,
                    dispersion_unit=dispersion_unit if dispersion_unit is
                                                       not None else
                                                       self.dispersion_unit)

    @property
    def dispersion(self):
        """Dispersion values."""
        if self._dispersion is None:
            _dispersion = np.arange(self.data.size)

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
                    _dispersion = np.logspace(crval, end, num)
                else:
                    _dispersion = np.arange(crval, end, cdelt)

                self._dispersion = np.ma.array(data=_dispersion,
                                               mask=self.mask)
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
        if type(other) == Layer or type(self) == Layer:
            this = self

            if isinstance(self, ModelLayer):
                this = other
                other = self

            # Make sure units are compatible
            if not other.unit.is_equivalent(
                    this.unit,
                    equivalencies=spectral_density(other.dispersion)):
                logging.error("Spectral data objects have incompatible units.")
                return

            # Check sampling; re-sample if incompatible
            other_samp = np.ones(this.data.data.value.shape,
                                 [('wlen', float), ('flux', float),
                                  ('ivar', float)])
            other_samp['wlen'] = other.dispersion.data.value
            other_samp['flux'] = other.data.data.value

            if isinstance(other, ModelLayer):
                cols = ('flux',)
            else:
                other_samp['ivar'] = other.uncertainty.data.value
                cols = ('flux', 'ivar')

            other_resamp = resample(other_samp, 'wlen',
                                      this.dispersion.data.value,
                                      cols)

            other_data = other._source.copy(
                other_resamp['flux'],
                uncertainty=other_resamp['ivar']
                            if not isinstance(other, ModelLayer) else None,
                mask=other._source.mask)
            other_data._dispersion = other_resamp['wlen']
        elif isinstance(self, ModelLayer) and isinstance(other, ModelLayer):
            # In the case where the model layers are from two separate data
            # layers, create a combine data layer
            # operator_char = dict(add='+', subtract='-', divide='-',
            #                      multiply='*')[operator]
            # final_layer = Layer.from_formula(
            #     "{} {} {}".format(self._parent.name, operator_char,
            #                       other._parent.name),
            #     [self._parent, other._parent])

            operator_func = dict(add='add', subtract='sub',
                                 divide='truediv', multiply='mult')[operator]
            final_model = getattr(self.model,
                                  "__{}__".format(operator_func))(other.model)

            final_names = []

            for smod in final_model._submodels:
                if smod._name in final_names:
                    split_name = re.findall(r"[^\W\d_]+|\d+", smod._name)
                    split_name[-1] = str(int(split_name[-1]) + 1)
                    smod._name = "".join(split_name)

                final_names.append(smod._name)

            final_mask = self.mask | other.mask

            result_model_layer = ModelLayer(model=final_model,
                                            source=self._parent._source,
                                            mask=~final_mask,
                                            parent=self._parent)

            return result_model_layer
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

        # Perform arithmetic operation
        operator = getattr(self._source, operator)

        # Check for WCS incompatibility
        if self._source.wcs != other_data.wcs:
            logging.warning("WCS objects are not equivalent; overriding "
                            "wcs information on 'other'.")
            tmp_wcs = other_data._wcs
            other_data._wcs = self._source.wcs
            result_source = operator(other_data, propagate_uncertainties=propagate)
            other_data._wcs = tmp_wcs
        else:
            result_source = operator(other_data, propagate_uncertainties=propagate)

        if result_source.dispersion_unit.is_unity():
            result_source._dispersion_unit = other_data.dispersion_unit

        # Create a layer from the source data object
        result_layer = Layer(result_source, parent=self._parent,
                             name=self.name)

        return result_layer

    @classmethod
    def from_formula(cls, formula, layers):
        if not formula:
            return

        layers = layers
        new_layer = cls._evaluate(layers, formula)

        if new_layer is None:
            logging.error(
                "Failed to create new layer from formula: {}".format(formula))
            return

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
        source_data = np.ma.array(
            Quantity(self._source.data,
                     unit=self._source.unit).to(self.units[1]),
            mask=self._source.mask)
        layer_data = np.ma.array(source_data, mask=~self._mask)

        return layer_data

    @property
    def unit(self):
        """Flux unit."""
        return self.data.data.unit

    @property
    def dispersion(self):
        """Dispersion quantity with mask applied."""
        source_dispersion = np.ma.array(
            Quantity(self._source.dispersion,
                     unit=self._source.dispersion_unit).to(self.units[0]),
            mask=self._source.mask)
        layer_dispersion = np.ma.array(source_dispersion, mask=~self._mask)

        return layer_dispersion

    @dispersion.setter
    def dispersion(self, value, unit=""):
        self._source._dispersion = value

    @property
    def uncertainty(self):
        """Flux uncertainty with mask applied."""
        source_uncertainty = np.ma.array(
            Quantity(self._source.uncertainty.array,
                     unit=self._source.unit).to(self.units[1]),
            mask=self._source.mask)
        layer_uncertainty = np.ma.array(source_uncertainty, mask=~self._mask)

        return layer_uncertainty

    @property
    def mask(self):
        """Mask for spectrum data."""
        return self._source.mask.astype(bool) & ~self._mask.astype(bool)

    @property
    def wcs(self):
        """WCS for spectrum data."""
        return self._source.wcs

    @property
    def meta(self):
        """Spectrum metadata."""
        return self._source.meta

    def convert(self, data, unit, equivalencies=None):
        mask = None

        if isinstance(data, np.ma.MaskedArray):
            mask = data.mask
            data = data.data.to(unit, equivalencies=equivalencies)
        elif isinstance(data, Quantity):
            data = data.to(unit, equivalencies=equivalencies)

        new_data = np.ma.array(data,
                               mask=mask if mask is not None else self.mask)

        return new_data


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
        self._data = self._model(self.dispersion.data.value)
        layer_data = np.ma.array(
            Quantity(self._data,
                     unit=self._source.unit).to(self.units[1]),
            mask=self._source.mask)

        return layer_data

    @property
    def roi_data(self):
        """
        Apply an ROI mask to the model layer data.
        """
        return np.ma.array(self.data, mask=~self._mask)

    @property
    def dispersion(self):
        """Dispersion quantity with mask applied."""
        layer_dispersion = np.ma.array(
            Quantity(self._source.dispersion,
                     unit=self._source.dispersion_unit).to(self.units[0]),
            mask=self._source.mask)

        return layer_dispersion

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
        self._data = self._model(self.dispersion.data.value)

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