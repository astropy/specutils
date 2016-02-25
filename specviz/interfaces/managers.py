from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# LOCAL
from .factories import DataFactory, ModelFactory, PlotFactory, FitterFactory
from .initializers import initialize
from ..core.comms import EventNode
from ..analysis import modeling
from ..third_party.py_expression_eval import Parser
from ..core.comms import Dispatch, DispatchHandle
from ..analysis.modeling import apply_model

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

        Dispatch.on_added_data.emit(data)

    def remove(self, data):
        self._members.remove(data)

        Dispatch.on_removed_data.emit(data)


class WindowManager(Manager):
    """
    Manages the association of layer objects with sub windows.
    """
    def __init__(self):
        super(WindowManager, self).__init__()
        self._members = {}

        DispatchHandle.setup(self)

    def add(self, layer, window):
        if window not in self._members:
            self._members[window] = [layer]
        else:
            self._members[window].append(layer)

        Dispatch.on_added_to_window.emit(layer=layer, window=window)

    def remove(self, layer, window=None):
        if window in self._members:
            if layer in self._members[window]:
                self._members[window].remove(layer)
        else:
            for window in self._members:
                if layer in self._members[window]:
                    self._members[window].remove(layer)
                    break

        Dispatch.on_removed_from_window.emit(layer=layer, window=window)

    def get(self, layer):
        for window in self._members:
            for l in self._members[window]:
                if layer == l:
                    return window

    def get_layers(self, window):
        return self._members.get(window, [])

class LayerManager(Manager):
    """
    Manages a set of layer objects.
    """
    def __init__(self):
        super(LayerManager, self).__init__()

    def new(self, data, mask=None, parent=None, name='', model=None):
        logging.info("Creating new layer: {}".format(name))

        new_layer = DataFactory.create_layer(data, mask, parent, name, model)
        self.add(new_layer)

        return new_layer

    def add(self, layer):
        self._members.append(layer)

        # Emit creation event
        Dispatch.on_added_layer.emit(layer)

    def remove(self, layer):
        """
        Remove specified layer from the manager.

        Parameters
        ----------
        layer : specviz.core.data.Layer
            Layer object to be removed.
        """
        self._members.remove(layer)

        # Emit removal event
        Dispatch.on_removed_layer.emit(layer)

    def copy(self, layer):
        new_layer = DataFactory.create_layer(layer._source)
        new_layer.dispersion = layer.dispersion

        return new_layer

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
        new_layer = self.new(new_data, parent=layer._parent)
        new_layer.dispersion = layer.dispersion.value

        return new_layer

    def update_model_layer(self, layer, model, fitter=None):
        new_model = modeling.apply_model(model, layer.dispersion,
                                         layer.data, fitter)
        # layer.set_model(new_model)

    def update_model_parameters(self, model, parameters):
        pass

    def add_from_formula(self, formula):
        if not formula:
            return

        new_layer = self._evaluate(self._members, formula)
        new_layer.name = "Resultant"
        new_layer._window = None
        self.add(new_layer)

        return new_layer

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

        for layer in layers:
            formula = formula.replace(layer.name, layer.name.replace(" ", "_"))

        expr = parser.parse(formula)

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

        result = parser.evaluate(expr.simplify({}).toString(),
                                 dict(pair for pair in zip(vars, sorted_layers)))

        return result


class ModelManager(Manager):
    """
    Manages a set of model objects.
    """
    def __init__(self):
        super(ModelManager, self).__init__()
        self._members = {}
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
        flux = data_layer.data[mask]
        dispersion = data_layer.dispersion[mask]
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

        if not formula:
            result = np.sum(models) if len(models) > 1 else models[0]
            return result

        return self._evaluate(models, formula)

    def transfer_models(self, old_layer, new_layer):
        mdls = self._members[old_layer]
        self._members[new_layer] = mdls

        self._members[old_layer] = []
        del self._members[old_layer]

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

        mask = parent_layer._mask
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


class PlotManager(Manager):
    """
    Manages all plots.
    """
    def __init__(self):
        super(PlotManager, self).__init__()
        self._members = {}

    def new(self, layer, window, unit=None, visible=False, style='line',
            pen=None):
        plot_container = PlotFactory.create_line_plot(layer, unit, visible,
                                                      style, pen)

        self.add(plot_container, layer, window)

        return plot_container

    def add(self, plot_container, layer, window):
        if window not in self._members:
            self._members[window] = [plot_container]
        else:
            self._members[window].append(plot_container)

        Dispatch.on_added_plot.emit(container=plot_container, window=window)

    def remove(self, layer, window=None):
        if window is None:
            for w in self._members:
                if layer in self._members:
                    window = w
                    break

        Dispatch.on_removed_plot.emit(layer=layer, window=window)

    def get_plots(self, window):
        return self._members[window]

    def get_plot_from_layer(self, layer, window):
        for container in self._members[window]:
            if container.layer == layer:
                return container

        logging.warning("No plots found for layer {} in window {}.".format(
            layer, window))

    def update_plots(self, container=None, layer=None):
        if container is not None:
            container.update()
        elif layer is not None:
            for container in [y for x in self._members for y in
                              self._members[x]]:
                if container._layer == layer:
                    container.update()

        Dispatch.on_updated_plot.emit(container=container, layer=layer)


data_manager = DataManager()
window_manager = WindowManager()
layer_manager = LayerManager()
model_manager = ModelManager()
plot_manager = PlotManager()
