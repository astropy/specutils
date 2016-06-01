from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# STDLIB
import logging

# LOCAL
from ..core.data import Data, Layer, ModelLayer
from ..core.plots import LinePlot
from ..analysis.models.spline import Spline1D
from ..analysis.models.blackbody import BlackBody

from ..ui.widgets.sub_windows import PlotSubWindow

# THIRD-PARTY
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
        new_data = Data.read(path, filter)

        return new_data

    @staticmethod
    def from_array(array):
        new_data = Data(array)

        return new_data

    @staticmethod
    def create_layer(data, mask=None, parent=None, name='', model=None):
        if not isinstance(data, Data):
            raise ValueError('Invalid Data object')

        mask = mask if mask is not None else np.ones(data.data.shape,
                                                     dtype=bool)

        if model is None:
            new_layer = Layer(data, mask, parent, name)
        else:
            new_layer = ModelLayer(model, data, mask, parent, name)

        return new_layer


#TODO  a base class for Model and Fitter classes might be of help here.

class ModelFactory(Factory):

    # Ideally we should be getting these classes from astropy directly and
    # transparently, instead of explicitly naming them here. This is basically
    # a maintenance issue: at each new release of astropy we should check if
    # new models became available, and existing models got deprecated.
    #
    # This might not be possible unless astropy itself somehow manages to
    # further subclass its Fittable1DModel class into spectrum-specific and
    # other types. Right now we have a mix of spectral models, galaxy surface
    # brightness models, and others, all lumped into a single type. Thus we
    # have to keep picking the spectral relevant types by hand for the time
    # being.
    all_models = {
        'Gaussian': models.Gaussian1D,
        'GaussianAbsorption': models.GaussianAbsorption1D,
        'Lorentz': models.Lorentz1D,
        'MexicanHat': models.MexicanHat1D,
        'Trapezoid': models.Trapezoid1D,
        'ExponentialCutoffPowerLaw': models.ExponentialCutoffPowerLaw1D,
        'BrokenPowerLaw': models.BrokenPowerLaw1D,
        'LogParabola': models.LogParabola1D,
        'PowerLaw': models.PowerLaw1D,
        'Linear': models.Linear1D,
        'Const': models.Const1D,
        'Redshift': models.Redshift,
        'Scale': models.Scale,
        'Shift': models.Shift,
        'Sine': models.Sine1D,
        'Voigt': models.Voigt1D,
        'Spline': Spline1D,
        'BlackBody': BlackBody,

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

        logging.error("No such model {}".format(name))


class FitterFactory(Factory):
    default_fitter = fitting.LevMarLSQFitter

    all_fitters = {
        'Levenberg-Marquardt': fitting.LevMarLSQFitter,
        'Simplex': fitting.SimplexLSQFitter,
        'SLSQP': fitting.SLSQPLSQFitter,
    }

    @classmethod
    def create_fitter(cls, name):
        name = str(name)

        if name in cls.all_fitters:
            return cls.all_fitters[name]

        logging.error("No such fitter {}".format(name))


class WindowFactory(Factory):
    """
    Produces sub window objects.
    """
    @classmethod
    def create_window(cls):
        sub_window = PlotSubWindow()

        return sub_window


class PlotFactory(Factory):
    """
    Produces plot container objects.
    """

    @classmethod
    def create_line_plot(cls, layer, unit=None, visible=False, style='line',
                         pen=None, err_pen=None):
        # User attempts to plot before loading file, do nothing
        if not isinstance(layer, Layer):
            return

        plot_data_item = pg.PlotDataItem(layer.dispersion.value,
                                         layer.data.value)

        plot_container = LinePlot(layer=layer, plot=plot_data_item,
                                  visible=visible, style=style,
                                  pen=pen, err_pen=err_pen)

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
