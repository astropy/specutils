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
