from __future__ import absolute_import, division, print_function

from ..core.data import Data, Layer
import numpy as np


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

    @staticmethod
    def create_layer(data, mask=None, parent=None):
        mask = mask or np.ones(data.data.shape, dtype=bool)
        new_layer = Layer(data, mask, parent)
        return new_layer
