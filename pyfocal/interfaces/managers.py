from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# LOCAL
from .factories import DataFactory, ModelFactory, PlotFactory, FitterFactory
from ..core.events import EventHook
from ..analysis import modeling
from ..third_party.py_expression_eval import Parser

# STDLIB
import logging

# THIRD-PARTY
import numpy as np


class Manager(object):
    """
    Manages sets of objects.
    """
    def __init__(self):
        self._members = []

        self.on_add = EventHook()
        self.on_remove = EventHook()


class DataManager(Manager):
    """
    Manages a set of data objects.
    """
    def __init__(self):
        super(DataManager, self).__init__()

    def load(self, path, filter):
        new_data = DataFactory.from_file(path, filter)
        self.add(new_data)

        return new_data

    def add(self, data):
        self._members.append(data)

    def remove(self, data):
        self._members.remove(data)


class LayerManager(Manager):
    """
    Manages a set of layer objects.
    """
    def __init__(self):
        super(LayerManager, self).__init__()

    def new_layer(self, data, mask=None, parent=None, window=None, name=''):
        logging.info("Creating new layer: {}".format(name))
        new_layer = DataFactory.create_layer(data, mask, parent, window, name)

        self._members.append(new_layer)

        # Emit creation event
        self.on_add.emit(new_layer)

        return new_layer

    def new_model_layer(self, model, data, mask, parent=None, window=None,
                        name=''):
        logging.info("Creating new model layer: {}".format(name))
        model_layer = DataFactory.create_model_layer(model, data, mask,
                                                     parent=parent,
                                                     window=window,
                                                     name=name)

        self._members.append(model_layer)

        # Emit creation event
        self.on_add.emit(model_layer)

        return model_layer

    def add_layer(self, layer):
        self._members.append(layer)

        self.on_add.emit(layer)

    def remove(self, layer):
        """
        Remove specified layer from the manager.

        Parameters
        ----------
        layer : pyfocal.core.data.Layer
            Layer object to be removed.
        """
        self._members.remove(layer)

        # Emit removal event
        self.on_remove.emit(layer)

    def get_window_layers(self, window):
        """
        Retrieve all children of the `SubWindow` object.
        """
        return [x for x in self._members if x._window == window]

    def add_from_model(self, layer, model, fitter=None, mask=None):
        if layer is None:
            logging.error("No layer selected from which to create a model.")
            return

        new_model = modeling.apply_model(model, layer.dispersion,
                                         layer.data, fitter)

        new_data = DataFactory.from_array(new_model(layer.dispersion.value))
        new_layer = self.new_layer(new_data,
                                   parent=layer._parent, window=layer._window)
        new_layer.dispersion = layer.dispersion.value

        return new_layer

    def update_model_layer(self, layer, model, fitter=None):
        new_model = modeling.apply_model(model, layer.dispersion,
                                         layer.data, fitter)
        # layer.set_model(new_model)

    def update_model_parameters(self, model, parameters):
        pass

    def _evaluate(self, layers, formula):
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
        expr = parser.parse(formula)

        # Extract variables
        vars = expr.variables()

        # List the models in the same order as the variables
        sorted_layers = [l for v in vars for l in layers if l.name == v]

        if len(sorted_layers) != len(vars):
            logging.error("Incorrect layer arithmetic formula: the number "
                          "of layers does not match the number of variables.")

        result = parser.evaluate(expr.simplify({}).toString(),
                                 dict(pair for pair in zip(vars, sorted_layers)))

        return result


class ModelManager(Manager):
    """
    Manages a set of model objects.
    """
    def __init__(self):
        super(ModelManager, self).__init__()
        # Model layer manager specified events
        self.on_add_model = EventHook()
        self.on_remove_model = EventHook()

        self._members = {}
        self.all_models = sorted(ModelFactory.all_models.keys())
        self.all_fitters = sorted(FitterFactory.all_fitters.keys())

    def new_model(self, layer, model_name):
        model = ModelFactory.create_model(model_name)()

        self.add(model, layer)

        self.on_add_model.emit(model)

        return model

    def add(self, model, layer):
        if layer not in self._members:
            self._members[layer] = [model]
        else:
            self._members[layer].append(model)

        self.on_add.emit(layer, model)

        return model

    def remove(self, layer, index=None):
        if index is not None:
            model_layer = self._members[layer].pop(index)
        else:
            del self._members[layer]

        self.on_remove.emit(layer)

    def _evaluate(self, models, formula):
        parser = Parser()
        expr = parser.parse(formula)

        # Extract variables
        vars = expr.variables()

        # List the models in the same order as the variables
        sorted_models = [m for v in vars for m in models if m.name == v]

        if len(sorted_models) != len(vars):
            logging.error("Incorrect model arithmetic formula: the number "
                          "of models does not match the number of variables.")

        result = parser.evaluate(expr.simplify({}).toString(),
                                 dict(pair for pair in zip(vars, sorted_models)))

        return result

    def get_model_layers(self, layer, no_keys=False):
        """
        Returns the `ModelLayer` objects associated with `Layer`.

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
            model_layers = []

            for k, v in self._members.items():
                if k._source == layer:
                    model_layers.append(*v)

            result += model_layers

        return result

    def get_compound_model(self, model_dict, formula=''):
        models = []

        for model in model_dict.keys():
            for i, param_name in enumerate(model.param_names):
                setattr(model, param_name, model_dict[model][i])

            models.append(model)

        if not formula:
            result = np.sum(models) if len(models) > 1 else models[0]
            return result

        return self._evaluate(models, formula)

    def transfer_models(self, old_layer, new_layer):
        mdls = self._members[old_layer]
        self._members[new_layer] = mdls

        self._members[old_layer] = []
        del self._members[old_layer]

    def update_model(self, model_layer, model_dict, formula=''):
        logging.info("ModelManager.update_model: {}".format(model_dict))
        model = self.get_compound_model(model_dict, formula)
        model_layer.model = model


class PlotManager(Manager):
    """
    Manages all plots.
    """
    def __init__(self):
        super(PlotManager, self).__init__()
        self._members = {}

    def new_line_plot(self, layer, parent, unit=None, visible=False,
                      style='line', pen=None):
        plot_container = PlotFactory.create_line_plot(layer, unit, visible,
                                                      style, pen)

        if parent not in self._members:
            self._members[parent] = [plot_container]
        else:
            self._members[parent].append(plot_container)

        return plot_container

    def get_plots(self, parent):
        return self._members[parent]

    def update_plots(self, window, layer):
        for container in self._members[window]:
            if container.layer == layer:
                container.update()


data_manager = DataManager()
layer_manager = LayerManager()
model_manager = ModelManager()
plot_manager = PlotManager()
