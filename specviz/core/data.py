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
from astropy.units import Unit, Quantity, spectral_density, spectral
from ..third_party.py_expression_eval import Parser
from specutils.core.generic import GenericSpectrum1D


class GenericSpectrum1DLayer(GenericSpectrum1D):
    """Class to handle layers in SpecViz."""
    def __init__(self, parent=None, layer_mask=None, *args, **kwargs):
        super(GenericSpectrum1DLayer, self).__init__(*args, **kwargs)
        self._parent = parent
        self._layer_mask = layer_mask

    @classmethod
    def from_parent(cls, parent, layer_mask=None):
        return cls(name=parent.name + " Layer", data=parent.data,
                   unit=parent.unit, uncertainty=parent.uncertainty,
                   mask=parent.mask, wcs=parent.wcs,
                   dispersion=parent.dispersion,
                   dispersion_unit=parent.dispersion_unit,
                   layer_mask=layer_mask, parent=parent,
                   copy=False)

    def from_self(self, name="", layer_mask=None):
        gen_spec = super(GenericSpectrum1DLayer, self).from_self(name=name)

        return self.from_parent(parent=gen_spec, layer_mask=layer_mask)

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
        """Flux quantity with mask applied. Returns a masked array
        containing a Quantity object."""
        data = np.ma.array(
            Quantity(self._data, unit=self.unit),
            mask=self.full_mask)

        return data

    @property
    def dispersion(self):
        """Dispersion quantity with mask applied. Returns a masked array
        containing a Quantity object."""
        dispersion = np.ma.array(
            Quantity(self._dispersion, unit=self.dispersion_unit),
            mask=self.full_mask)

        return dispersion

    @GenericSpectrum1D.uncertainty.getter
    def uncertainty(self):
        """Flux uncertainty with mask applied. Returns a masked array
        containing a Quantity object."""
        uncertainty = np.ma.array(
            Quantity(self._uncertainty.array, unit=self.unit),
            mask=self.full_mask)

        return uncertainty

    @property
    def layer_mask(self):
        """Mask applied from an ROI."""
        if self._layer_mask is None:
            self._layer_mask = np.ones(self._data.shape).astype(bool)

        return self._layer_mask

    @property
    def full_mask(self):
        """Mask for spectrum data."""
        return self.mask.astype(bool) | ~self.layer_mask.astype(bool)

    def set_units(self, disp_unit, data_unit):
        if self.dispersion_unit.is_equivalent(disp_unit,
                                              equivalencies=spectral()):
            self._dispersion = self.dispersion.data.to(
                disp_unit, equivalencies=spectral()).value

            # Finally, change the unit
            self.dispersion_unit = disp_unit
        else:
            logging.warning("Units are not compatible.")

        if self.unit.is_equivalent(data_unit,
                                   equivalencies=spectral_density(
                                       self.dispersion.data)):
            self._data = self.data.data.to(
                data_unit, equivalencies=spectral_density(
                    self.dispersion.data)).value

            self._uncertainty = self._uncertainty.__class__(
                self.uncertainty.data.to(
                    data_unit, equivalencies=spectral_density(
                        self.dispersion.data)).value)

            # Finally, change the unit
            self._unit = data_unit
        else:
            logging.warning("Units are not compatible.")


class GenericSpectrum1DModelLayer(GenericSpectrum1DLayer):
    """A layer for spectrum with a model applied."""
    def __init__(self, model=None, *args, **kwargs):
        super(GenericSpectrum1DModelLayer, self).__init__(*args, **kwargs)
        self._model = model

    @classmethod
    def from_parent(cls, parent, model=None, layer_mask=None):
        return cls(name=parent._parent.name + " Model Layer",
                   data=model(parent.dispersion.data.value),
                   unit=parent.unit,
                   uncertainty=parent.uncertainty.__class__(np.zeros(
                       parent.data.shape)),
                   mask=parent.mask, wcs=parent.wcs,
                   dispersion=parent.dispersion,
                   dispersion_unit=parent.dispersion_unit,
                   layer_mask=layer_mask, parent=parent, model=model,
                   copy=False)

    @property
    def model(self):
        """Spectrum model."""
        return self._model

    @model.setter
    def model(self, value):
        self._model = value
        self._data = self._model(self.dispersion.data.value)

    @classmethod
    def from_formula(cls, models, formula):
        result_model = cls._evaluate(models, formula)

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