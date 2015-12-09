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


data_manager = DataManager()