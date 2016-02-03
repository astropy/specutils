from .factories import DataFactory, ModelFactory, PlotFactory
from ..core.events import EventHook
from ..analysis import modeling

from py_expression_eval import Parser
import logging
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

    def new(self, data, mask=None, parent=None):
        new_layer = DataFactory.create_layer(data, mask, parent)

        return new_layer

    def add(self, data, mask=None, parent=None):
        new_layer = self.new(data, mask, parent)
        self._members.append(new_layer)

        # Emit creation event
        self.on_add.emit(new_layer)

        return new_layer

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

    def get_sub_window_layers(self, sub_window):
        """
        Retrieve all children of the `SubWindow` object.
        """
        return [x for x in self._members if x._parent == sub_window]

    def add_from_model(self, layer, model, fitter=None, mask=None):
        if layer is None:
            logging.error("No layer selected from which to create a model.")
            return

        new_model = modeling.apply_model(model, layer.dispersion,
                                         layer.data, fitter)

        new_data = DataFactory.from_array(new_model(layer.dispersion.value))
        new_layer = self.add(new_data, parent=layer._parent)
        new_layer.dispersion = layer.dispersion.value
        # new_layer.set_model(new_model)

        return new_layer

    def update_model_layer(self, layer, model, fitting=None):
        new_model = modeling.apply_model(model)
        layer.set_model(new_model)

    def update_model_parameters(self, model, parameters):
        pass


class ModelManager(Manager):
    """
    Manages a set of model objects.
    """
    def __init__(self):
        super(ModelManager, self).__init__()
        self._members = {}
        self.all_models = ModelFactory.all_models.keys()

    def add(self, layer, name):
        model = ModelFactory.create_model(name)()

        if layer not in self._members:
            self._members[layer] = [model]
        else:
            self._members[layer].append(model)

        self.on_add.emit(layer, model, name)

        return model

    def remove(self, layer, index):
        model = self._members[layer].pop(index)

        self.on_remove.emit(layer, model)

    def _evaluate(self, models, formula):
        parser = Parser()
        expr = parser.parse(formula)
        vars = expr.variables()
        result = parser.evaluate(expr.simplify({}).toString(),
                                 dict(pair for pair in zip(vars, models)))

        return result

    def get_layer_models(self, layer):
        return self._members.get(layer, [])

    def get_compound_model(self, layer, model_dict, formula=''):
        models = []

        for model in self._members[layer]:
            for i, param_name in enumerate(model.param_names):
                setattr(model, param_name, model_dict[model][i])

            models.append(model)

        if not formula:
            print(models)
            result = np.sum(models) if len(models) > 1 else models[0]
            return result

        return self._evaluate(models, formula)


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



data_manager = DataManager()
layer_manager = LayerManager()
model_manager = ModelManager()
plot_manager = PlotManager()
