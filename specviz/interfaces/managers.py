from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# LOCAL
from .factories import ModelFactory, FitterFactory
from .initializers import initialize
from ..analysis import modeling
from ..third_party.py_expression_eval import Parser
from ..core.comms import Dispatch, DispatchHandle
from ..analysis.modeling import apply_model

# STDLIB
import logging
from collections import OrderedDict

# THIRD-PARTY
import numpy as np


class Manager(object):
    """
    Manages sets of objects.
    """
    def __init__(self):
        self._members = []


class ModelManager(Manager):
    """
    Manages a set of model objects.
    """
    def __init__(self):
        super(ModelManager, self).__init__()
        self._members = OrderedDict()
        self.all_models = sorted(ModelFactory.all_models.keys())
        self.all_fitters = sorted(FitterFactory.all_fitters.keys())

    def new(self, model_name, layer, mask):
        model = ModelFactory.create_model(model_name)()

        # if a model layer is selected, get data from its parent layer.
        data_layer = layer

        if hasattr(layer, '_model'):
            data_layer = layer._parent

        # initialize model with sensible parameter values.
        # mask must be provided by caller, since ROIs may
        # have been re-defined on the plot since last time
        # the layer was updated.
        flux = data_layer.data[mask[data_layer._mask]]
        dispersion = data_layer.dispersion[mask[data_layer._mask]]

        initialize(model, dispersion, flux)

        self.add(model, layer)

        return model

    def add(self, model, layer):
        if layer not in self._members:
            self._members[layer] = [model]
        else:
            self._members[layer].append(model)

        # Emit event
        Dispatch.on_added_model.emit(model, layer)

        return model

    def remove(self, layer, model=None, index=None):
        if layer not in self._members:
            return

        if index is not None:
            model = self._members[layer].pop(index)
        elif model is not None:
            if model in self._members[layer]:
                self._members[layer].remove(model)
        else:
            del self._members[layer]

        Dispatch.on_removed_model.emit(model=model, layer=layer)

    @classmethod
    def evaluate(self, models, formula):
        parser = Parser()
        expr = parser.parse(formula)

        # Extract variables
        vars = expr.variables()

        # List the models in the same order as the variables
        sorted_models = [m for v in vars for m in models if m.name == v]

        if len(sorted_models) != len(vars):
            logging.error("Incorrect model arithmetic formula: the number "
                          "of models does not match the number of variables.")
            return

        result = parser.evaluate(expr.simplify({}).toString(),
                                 dict(pair for pair in zip(vars, sorted_models)))

        return result

    def get_models(self, layer, no_keys=False):
        """
        Returns the astropy model objects associated with `Layer`.

        Parameters
        ----------
        layer : core.data.Layer
            The `Layer` object the `ModelLayer`s are associated with.
        no_keys : bool, optional
            Whether to exclude top-level `ModelLayer` objects (`True`).

        Returns
        -------
        result : list
            List of `ModelLayer` objects.
        """
        result = self._members.get(layer, [])

        if not no_keys:
            models = []

            for k, v in self._members.items():
                if k._source == layer:
                    models.append(*v)

            result += models

        return result

    def get_compound_model(self, model_dict, formula=''):
        models = []

        for model in model_dict:
            for i, param_name in enumerate(model.param_names):
                setattr(model, param_name, model_dict[model][i])

            models.append(model)

        if formula:
            result = self.evaluate(models, formula)

            if result is not None:
                return result

        return np.sum(models) if len(models) > 1 else models[0]

    def transfer_models(self, old_layer, new_layer):
        mdls = self._members[old_layer]
        self._members[new_layer] = mdls

        self._members[old_layer] = []
        del self._members[old_layer]

        for model in self._members[new_layer]:
            Dispatch.on_added_model.emit(model, new_layer)

    def update_model(self, layer, model_inputs, formula='', mask=None):
        model = self.get_compound_model(model_inputs, formula)

        if hasattr(layer, '_model'):
            layer.model = model

        if mask is not None:
            layer._mask = mask

        Dispatch.on_updated_model.emit(model=model)

    def update_model_parameters(self, model, model_inputs):
        for model in model_inputs:
            model.parameters = model_inputs[model]

    def fit_model(self, layer, fitter_name):
        if not hasattr(layer, 'model'):
            logging.warning("This layer has no model to fit.")
            return

        # When fitting, the selected layer is a ModelLayer, thus
        # the data to be fitted resides in the parent
        parent_layer = layer._parent

        if parent_layer is None:
            return

        # While the data comes from the parent, the mask from the model
        # layer is the actual data that needs to be fit
        mask = layer._mask
        flux = parent_layer.data[mask]
        dispersion = parent_layer.dispersion[mask]
        model = layer.model

        # If the number of parameters is greater than the number of data
        # points, bail
        if len(model.parameters) > flux.size:
            logging.warning("Unable to perform fit; number of parameters is "
                            "greater than the number of data points.")
            return

        fitted_model = apply_model(model, dispersion, flux,
                                   fitter_name=fitter_name)

        # Update original model with new values from fitted model
        if hasattr(fitted_model, '_submodels'):
            for i in range(len(fitted_model._submodels)):
                for pname in model._submodels[i].param_names:
                    value = getattr(fitted_model, "{}_{}".format(pname, i))
                    setattr(model._submodels[i], pname, value.value)
                    setattr(model[i], pname, value.value)
        else:
            for pname in model.param_names:
                value = getattr(fitted_model, "{}".format(pname))
                setattr(model, pname, value.value)


        # update GUI with fit results
        Dispatch.on_updated_model.emit(model=model)

        return layer

model_manager = ModelManager()