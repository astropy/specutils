from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ..core.data import Data, Layer, ModelLayer
from ..core.containers import PlotContainer

import pyqtgraph as pg
import numpy as np
from astropy.modeling import models, fitting


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
        print("Creating from file.")
        new_data = Data.read(path, filter)
        return new_data

    @staticmethod
    def from_array(array):
        new_data = Data(array)
        return new_data

    @staticmethod
    def create_layer(data, mask=None, parent=None, window=None):
        mask = mask if mask is not None else np.ones(data.data.shape,
                                                     dtype=bool)
        new_layer = Layer(data, mask, parent, window)
        print("DataFactory.create_layer: {}, {}".format(new_layer._parent,
                                                     parent))
        return new_layer

    @staticmethod
    def create_model_layer(layer, model, parent=None):
        new_model_layer = ModelLayer(layer, model, parent)

        return new_model_layer


#TODO  a base class for Model and Fitter classes might be of help here.

class ModelFactory(Factory):

    # Any reason for just these three?
    # all_models = {
    #     "Gaussian1D": models.Gaussian1D,
    #     "Linear1D": models.Linear1D,
    #     "Const1D": models.Const1D
    # }

    # Ideally we should be getting these classes from astropy directly and
    # transparently, instead of explicitly naming them here. This is basically
    # a maintenance issue: at each new release of astropy we should check if
    # new models became available, and existing models got deprecated.
    all_models = {
        'Gaussian1D': models.Gaussian1D,
        'GaussianAbsorption1D': models.GaussianAbsorption1D,
        'Lorentz1D': models.Lorentz1D,
        'MexicanHat1D': models.MexicanHat1D,
        'Trapezoid1D': models.Trapezoid1D,
        'ExponentialCutoffPowerLaw1D': models.ExponentialCutoffPowerLaw1D,
        'BrokenPowerLaw1D': models.BrokenPowerLaw1D,
        'LogParabola1D': models.LogParabola1D,
        'PowerLaw1D': models.PowerLaw1D,
        'Linear1D': models.Linear1D,
        'Const1D': models.Const1D,
        'Redshift': models.Redshift,
        'Scale': models.Scale,
        'Shift': models.Shift,
        'Sine1D': models.Sine1D,
        # polynomials have to be handled separately. Their calling sequence
        # is incompatible with the Fittable1DModel interface, and they run
        # under a linear minimization algorithm, as opposed to the non-linear
        # minimization used with Fittable1DModel types.
        # 'Chebyshev1D': models.Chebyshev1D,
        # 'Legendre1D': models.Legendre1D,
        # 'Polynomial1D': models.Polynomial1D,
    }

    @classmethod
    def create_model(cls, name):
        name = str(name)

        if name in cls.all_models:
            return cls.all_models[name]

        print("No such model {}".format(name))


class FitterFactory(Factory):
    default_fitter = fitting.LevMarLSQFitter

    all_fitters = {
        'Levenberg-Marquardt': fitting.LevMarLSQFitter,
    }

    @classmethod
    def create_fitter(cls, name):
        name = str(name)

        if name in cls.all_fitters:
            return cls.all_fitters[name]

        print("No such fitter {}".format(name))


class PlotFactory(Factory):
    """
    Produces plot container objects.
    """

    @classmethod
    def create_line_plot(cls, layer, unit=None, visible=False, style='line',
                         pen=None, err_pen=None):
        plot_container = PlotContainer(layer=layer, visible=visible,
                                       style=style, pen=pen, err_pen=err_pen)

        plot_data_item = pg.PlotDataItem(plot_container.layer.dispersion.value,
                                         plot_container.layer.data.value)

        plot_container.plot = plot_data_item

        if plot_container.layer.uncertainty is not None:
            # err_top = pg.PlotDataItem(
            #     plot_container.layer.dispersion.value,
            #     plot_container.layer.data.value +
            #     plot_container.layer.uncertainty.array * 0.5)
            # err_btm = pg.PlotDataItem(
            #     plot_container.layer.dispersion.value,
            #     plot_container.layer.data.value -
            #     plot_container.layer.uncertainty.array * 0.5)
            #
            # plot_error_item = pg.FillBetweenItem(err_top, err_btm, 'r')

            plot_error_item = pg.ErrorBarItem(
                x=plot_container.layer.dispersion.value,
                y=plot_container.layer.data.value,
                height=plot_container.layer.uncertainty.array,
            )
            plot_container.error = plot_error_item

        return plot_container
