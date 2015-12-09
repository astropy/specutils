from __future__ import absolute_import, division, print_function

from ..core.data import Data


class Factory(object):
    """
    Responsible for creation of objects.
    """


class DataFactory(Factory):
    """
    Produces data objects.
    """
    def __init__(self):
        pass

    @staticmethod
    def from_file(path, filter):
        new_data = Data.read(path, filter)
        return new_data
