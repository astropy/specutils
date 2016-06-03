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

    def transfer_models(self, old_layer, new_layer):
        mdls = self._members[old_layer]
        self._members[new_layer] = mdls

        self._members[old_layer] = []
        del self._members[old_layer]

        for model in self._members[new_layer]:
            Dispatch.on_added_model.emit(model, new_layer)

    def update_model_parameters(self, model, model_inputs):
        for model in model_inputs:
            model.parameters = model_inputs[model]

model_manager = ModelManager()