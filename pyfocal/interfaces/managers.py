from __future__ import absolute_import, division, print_function

from .factories import DataFactory


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

    def get_sub_window_layers(self, sub_window):
        """
        Retrieve all children of the `SubWindow` object.
        """
        return [x for x in self._members if x._parent == sub_window]


data_manager = DataManager()
layer_manager = LayerManager()
