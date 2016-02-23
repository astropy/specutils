# This module is used to initialize spectral models to the
# data at hand. This is used by model-fitting code that has
# to create spectral model instances with sensible parameter
# values such that they can be used as first guesses by the
# fitting algorithms.

import numpy as np


def _get_model_name(model):
    class_string = str(model.__class__)
    return class_string.split('\'>')[0].split(".")[-1]


# Initialization that is specific to the Linear1D model.
class _Linear1DInitializer(object):

    def initialize(self, instance, x, y):

        y_range = np.max(y) - np.min(y)
        x_range = x[len(x) - 1] - x[0]
        slope = y_range / x_range
        y0 = y[0]

        instance.slope.value = slope
        instance.intercept.value = y0

        return instance


# Initialization that is specific to the Const1D model.
class _Const1DInitializer(object):

    def initialize(self, instance, x, y):
        average = np.mean(y)
        instance.amplitude.value = average
        return instance


# Initialization that is applicable to all "line profile"
# models, that is, models that have an amplitude, a width,
# and a defined position in wavelength space.
class _LineProfile1DInitializer(object):

    def __init__(self, factor=1.0):
        self._factor = factor

    def initialize(self, instance, x, y):

        y_range = np.max(y) - np.min(y)
        x_range = x[len(x) - 1] - x[0]
        position = x_range / 2.0 + x[0]
        width = x_range / 50.

        name = _get_model_name(instance)

        _setattr(instance, name, 'amplitude', y_range * self._factor)
        _setattr(instance, name, 'position', position)
        _setattr(instance, name, 'width', width)

        return instance


# Sets parameter value by mapping parameter name to model type.
# Prevents the parameter value setting to be stopped on its tracks
# by non-existent model names or parameter names.
def _setattr(instance, mname, pname, value):
    try:
        # conflicts happen when we try to set a parameter value
        # as a Quantity. Use the raw value instead.
        setattr(instance, _p_names[mname][pname], value.value)
    except KeyError:
        pass


# This associates each initializer to its corresponding spectral model.
# Some models are not really line profiles, but their parameter names
# and roles are the same as in a typical line profile, so they can be
# initialized in the same way.
_initializers = {
    'Beta1D':                     _LineProfile1DInitializer(),
    'Box1D':                      _LineProfile1DInitializer(),
    'Const1D':                    _Const1DInitializer(),
    'Gaussian1D':                 _LineProfile1DInitializer(),
    'GaussianAbsorption1D':       _LineProfile1DInitializer(),
    'Linear1D':                   _Linear1DInitializer(),
    'Lorentz1D':                  _LineProfile1DInitializer(),
    'Voigt1D':                    _LineProfile1DInitializer(),
    'MexicanHat1D':               _LineProfile1DInitializer(),
    'Trapezoid1D':                _LineProfile1DInitializer(),
    'PowerLaw1D':                 _LineProfile1DInitializer(factor=0.5),
    'BrokenPowerLaw1D':           _LineProfile1DInitializer(factor=0.5),
    'ExponentialCutoffPowerLaw1D':_LineProfile1DInitializer(factor=0.5),
    'LogParabola1D':              _LineProfile1DInitializer(factor=0.5),
}

# Models can have parameter names that are similar amongst them, but not quite the same.
_p_names = {
    'Gaussian1D':                 {'amplitude': 'amplitude', 'position': 'mean', 'width': 'stddev'},
    'GaussianAbsorption1D':       {'amplitude': 'amplitude', 'position': 'mean', 'width': 'stddev'},
    'Beta1D':                     {'amplitude': 'amplitude', 'position': 'x_0'},
    'Lorentz1D':                  {'amplitude': 'amplitude', 'position': 'x_0', 'width': 'fwhm'},
    'Voigt1D':                    {'amplitude': 'amplitude', 'position': 'x_0', 'width': 'fwhm_G'},
    'Box1D':                      {'amplitude': 'amplitude', 'position': 'x_0', 'width': 'width'},
    'MexicanHat1D':               {'amplitude': 'amplitude', 'position': 'x_0', 'width': 'sigma'},
    'Trapezoid1D':                {'amplitude': 'amplitude', 'position': 'x_0', 'width': 'width'},
    'PowerLaw1D':                 {'amplitude': 'amplitude', 'position': 'x_0'},
    'BrokenPowerLaw1D':           {'amplitude': 'amplitude', 'position': 'x_break'},
    'ExponentialCutoffPowerLaw1D':{'amplitude': 'amplitude', 'position': 'x_0'},
    'LogParabola1D':              {'amplitude': 'amplitude', 'position': 'x_0'},
    }


# Main entry point. X and Y are for now Quantity arrays with the
# independent and dependent variables. It's assumed X values
# are stored in increasing order in the array.
def initialize(instance, x, y):
    if x is None or y is None:
        return instance

    name = _get_model_name(instance)

    try:
        initializer = _initializers[name]

        return initializer.initialize(instance, x, y)

    except KeyError:
        return instance
