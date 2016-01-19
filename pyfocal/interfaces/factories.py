from ..core.data import Data, Layer
from ..core.containers import PlotContainer

import pyqtgraph as pg
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
    def from_array(array):
        new_data = Data(array)
        return new_data

    @staticmethod
    def create_layer(data, mask=None, parent=None):
        mask = mask if mask is not None else np.ones(data.data.shape,
                                                     dtype=bool)
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


class PlotFactory(Factory):
    """
    Produces plot container objects.
    """

    @classmethod
    def create_line_plot(cls, layer, unit=None, visible=False, style='line',
                         pen=None):
        plot_container = PlotContainer(layer=layer, visible=visible,
                                       style=style, pen=pen)

        plot_data_item = pg.PlotDataItem(plot_container.dispersion.value,
                                         plot_container.data.value)
        print("Creating plot")
        plot_container.plot = plot_data_item

        return plot_container