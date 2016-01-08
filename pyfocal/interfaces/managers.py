from .factories import DataFactory, ModelFactory
from ..core.events import EventHook

from py_expression_eval import Parser


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

    def remove(self, data):
        self._members.remove(data)


class LayerManager(Manager):
    """
    Manages a set of layer objects.
    """
    def __init__(self):
        super(LayerManager, self).__init__()

    def new(self, data, mask=None, sub_window=None):
        new_layer = DataFactory.create_layer(data, mask, sub_window)
        self._members.append(new_layer)
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

    def get_sub_window_layers(self, sub_window):
        """
        Retrieve all children of the `SubWindow` object.
        """
        return [x for x in self._members if x._parent == sub_window]


class ModelManager(Manager):
    """
    Manages a set of model objects.
    """
    on_add_model = EventHook()
    on_remove_model = EventHook()

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

        self.on_add_model.emit(layer, model)

        return model

    def remove(self, layer, index):
        model = self._members[layer].pop(index)

        self.on_remove_model.emit(layer, model)

    def evaluate(self, layer, formula):
        parser = Parser()
        expr = parser.parse(formula)
        vars = expr.variables()
        mdls = self._members[layer]
        result = parser.evaluate(expr.simplify({}).toString(),
                                 dict(pair for pair in zip(vars, mdls)))

        return result

    def get_layer_models(self, layer):
        return self._members.get(layer, [])


data_manager = DataManager()
layer_manager = LayerManager()
model_manager = ModelManager()
