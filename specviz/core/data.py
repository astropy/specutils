"""This module handles spectrum data objects."""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ..analysis.utils import resample
from .mixins import LayerArithmeticMixin

# STDLIB
import logging
logging.basicConfig(level=logging.INFO)
import re

# THIRD-PARTY
import numpy as np
from astropy.units import Unit, Quantity
from ..third_party.py_expression_eval import Parser
from specutils.core.generic import GenericSpectrum1D


class GenericSpectrum1DLayer(GenericSpectrum1D):
    """Class to handle layers in SpecViz."""
    def __init__(self, parent=None, *args, **kwargs):
        super(GenericSpectrum1DLayer, self).__init__(*args, **kwargs)
        self._parent = parent

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

    @property
    def data(self):
        """Flux quantity with mask applied.  Returns a masked array
        containing a Quantity object."""
        data = np.ma.array(
            Quantity(self._data, unit=self.unit),
            mask=self.full_mask)

        return data

    @property
    def dispersion(self):
        """Dispersion quantity with mask applied.  Returns a masked array
        containing a Quantity object."""
        dispersion = np.ma.array(
            Quantity(self._dispersion, unit=self.dispersion_unit),
            mask=self.full_mask)

        return dispersion

    @property
    def uncertainty(self):
        """Flux uncertainty with mask applied. Returns a masked array
        containing a Quantity object."""
        uncertainty = np.ma.array(
            Quantity(self._uncertainty, unit=self.unit),
            mask=self.full_mask)

        return uncertainty

    @property
    def full_mask(self):
        """Mask for spectrum data."""
        return self.mask.astype(bool) & ~self._layer_mask.astype(bool)

    def convert(self, data, unit, equivalencies=None):
        mask = None

        if isinstance(data, np.ma.MaskedArray):
            mask = data.mask
            data = data.data.to(unit, equivalencies=equivalencies)
        elif isinstance(data, Quantity):
            data = data.to(unit, equivalencies=equivalencies)

        new_data = np.ma.array(data,
                               mask=mask if mask is not None else self.full_mask)

        return new_data


class GenericSpectrum1DModelLayer(GenericSpectrum1DLayer):
    """A layer for spectrum with a model applied."""
    @property
    def roi_data(self):
        """
        Apply an ROI mask to the model layer data.
        """
        return np.ma.array(self.data, mask=~self._layer_mask)

    @property
    def uncertainty(self):
        """
        Models do not need to contain uncertainties; override parent
        class method.
        """
        return np.zeros(self.data.shape)

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


if __name__ == '__main__':
    flux1, disp1 = np.random.sample(100), np.arange(100)
    data1 = GenericSpectrum1D(flux1, dispersion=disp1)

    flux2, disp2 = np.random.sample(100), np.arange(100)
    data2 = GenericSpectrum1D(flux2, dispersion=disp2)

    layer1, layer2 = GenericSpectrum1DLayer(data1), GenericSpectrum1DLayer(data2)