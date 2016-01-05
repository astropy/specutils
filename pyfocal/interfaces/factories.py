from __future__ import absolute_import, division, print_function

from ..core.data import Data, Layer

import numpy as np
from astropy.modeling import models


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


class ModelFactory(Factory):
    all_models = {
        "Gaussian1D": models.Gaussian1D,
        "Linear1D": models.Linear1D,
        "Const1D": models.Const1D
    }

    @classmethod
    def create_model(cls, name):
        name = str(name)

        if name in cls.all_models:
            return cls.all_models[name]

        print("No such model {}".format(name))
