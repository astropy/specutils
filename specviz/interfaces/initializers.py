"""
This module is used to initialize spectral models to the
data at hand.

This is used by model-fitting code that has
to create spectral model instances with sensible parameter
values such that they can be used as first guesses by the
fitting algorithms.
"""
import numpy as np

from ..analysis.models import spline, blackbody

__all__ = [
    'initialize'
]

AMPLITUDE = 'amplitude'
POSITION  = 'position'
WIDTH     = 'width'


def _get_model_name(model):
    class_string = str(model.__class__)
    return class_string.split('\'>')[0].split(".")[-1]


class _Linear1DInitializer(object):
    """
    Initialization that is specific to the Linear1D model.

    Notes
    -----
    In a way, we need this specialized initializer because
    the linear 1D model is more like a kind of polynomial.
    It doesn't mesh well with other non-linear models.
    """
    def initialize(self, instance, x, y):
        """
        Initialize the model

        Parameters
        ----------
        instance: `~astropy.modeling.models`
            The model to initialize.

        x, y: numpy.ndarray
            The data to use to initialize from.

        Returns
        -------
        instance: `~astropy.modeling.models`
            The initialized model.
        """

        # y_range = np.max(y) - np.min(y)
        # x_range = x[-1] - x[0]
        # slope = y_range / x_range
        # y0 = y[0]

        y_mean = np.mean(y)

        instance.slope.value = 0.0
        instance.intercept.value = y_mean

        return instance


class _WideBand1DInitializer(object):
    """
    Initialization that is applicable to all "wide band"
    models

    A "wide band" model is one that has an amplitude and
    a position in wavelength space where this amplitude
    is defined.

    Parameters
    ----------
    factor: float
        The scale factor to apply to the amplitutde
    """
    def __init__(self, factor=1.0):
        self._factor = factor

    def initialize(self, instance, x, y):
        """
        Initialize the model

        Parameters
        ----------
        instance: `~astropy.modeling.models`
            The model to initialize.

        x, y: numpy.ndarray
            The data to use to initialize from.

        Returns
        -------
        instance: `~astropy.modeling.models`
            The initialized model.
        """
        y_mean = np.mean(y)
        x_range = x[-1] - x[0]
        position = x_range / 2.0 + x[0]

        name = _get_model_name(instance)

        _setattr(instance, name, AMPLITUDE, y_mean * self._factor)
        _setattr(instance, name, POSITION, position)

        return instance


class _LineProfile1DInitializer(object):
    """
    Initialization that is applicable to all "line profile"
    models.

    A "line profile" model is one that has an amplitude, a width,
    and a defined position in wavelength space.

    Parameters
    ----------
    factor: float
        The scale factor to apply to the amplitutde
    """
    def __init__(self, factor=1.0):
        self._factor = factor

    def initialize(self, instance, x, y):
        """
        Initialize the model

        Parameters
        ----------
        instance: `~astropy.modeling.models`
            The model to initialize.

        x, y: numpy.ndarray
            The data to use to initialize from.

        Returns
        -------
        instance: `~astropy.modeling.models`
            The initialized model.
        """

        # X centroid estimates the position
        centroid = np.sum(x * y) / np.sum(y)

        # width can be estimated by the weighted
        # 2nd moment of the X coordinate.
        dx = x - np.mean(x)
        width = np.sqrt(np.sum((dx * dx) * y) / np.sum(y))

        # amplitude is derived from area.
        sum = np.sum(y - np.min(y))
        height = sum / (width * np.sqrt(2 * np.pi) )

        name = _get_model_name(instance)

        _setattr(instance, name, AMPLITUDE, height * self._factor)
        _setattr(instance, name, POSITION, centroid)
        _setattr(instance, name, WIDTH, width)

        return instance


def _setattr(instance, mname, pname, value):
    """
    Sets parameter value by mapping parameter name to model type.

    Prevents the parameter value setting to be stopped on its tracks
    by non-existent model names or parameter names.

    Parameters
    ----------
    instance: `~astropy.modeling.models`
        The model to initialize.

    mname: str
        Model name.

    pname: str
        Parameter name.

    value: any
        The value to assign.
    """
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
    'Beta1D':                     _WideBand1DInitializer,
    'Const1D':                    _WideBand1DInitializer,
    'PowerLaw1D':                 _WideBand1DInitializer,
    'BrokenPowerLaw1D':           _WideBand1DInitializer,
    'ExponentialCutoffPowerLaw1D':_WideBand1DInitializer,
    'LogParabola1D':              _WideBand1DInitializer,
    'Box1D':                      _LineProfile1DInitializer,
    'Gaussian1D':                 _LineProfile1DInitializer,
    'GaussianAbsorption1D':       _LineProfile1DInitializer,
    'Lorentz1D':                  _LineProfile1DInitializer,
    'Voigt1D':                    _LineProfile1DInitializer,
    'MexicanHat1D':               _LineProfile1DInitializer,
    'Trapezoid1D':                _LineProfile1DInitializer,
    'Linear1D':                   _Linear1DInitializer,
    'Spline1D':                   spline.Spline1DInitializer,
    'BlackBody':                  blackbody.BlackBodyInitializer,
}

# Models can have parameter names that are similar amongst them, but not quite the same.
# This maps the standard names used in the code to the actual names used by astropy.
_p_names = {
    'Gaussian1D':                 {AMPLITUDE:'amplitude',  POSITION:'mean', WIDTH:'stddev'},
    'GaussianAbsorption1D':       {AMPLITUDE:'amplitude',  POSITION:'mean', WIDTH:'stddev'},
    'Lorentz1D':                  {AMPLITUDE:'amplitude',  POSITION:'x_0',  WIDTH:'fwhm'},
    'Voigt1D':                    {AMPLITUDE:'amplitude_L',POSITION:'x_0',  WIDTH:'fwhm_G'},
    'Box1D':                      {AMPLITUDE:'amplitude',  POSITION:'x_0',  WIDTH:'width'},
    'MexicanHat1D':               {AMPLITUDE:'amplitude',  POSITION:'x_0',  WIDTH:'sigma'},
    'Trapezoid1D':                {AMPLITUDE:'amplitude',  POSITION:'x_0',  WIDTH:'width'},
    'Beta1D':                     {AMPLITUDE:'amplitude',  POSITION:'x_0'},
    'PowerLaw1D':                 {AMPLITUDE:'amplitude',  POSITION:'x_0'},
    'ExponentialCutoffPowerLaw1D':{AMPLITUDE:'amplitude',  POSITION:'x_0'},
    'LogParabola1D':              {AMPLITUDE:'amplitude',  POSITION:'x_0'},
    'BrokenPowerLaw1D':           {AMPLITUDE:'amplitude',  POSITION:'x_break'},
    'Const1D':                    {AMPLITUDE:'amplitude'},
    }


def initialize(instance, x, y):
    """
    Initialize given model.

    X and Y are for now Quantity arrays with the
    independent and dependent variables. It's assumed X values
    are stored in increasing order in the array.

    Parameters
    ----------
    instance: `~astropy.modeling.models`
        The model to initialize.

    x, y: numpy.ndarray
        The data to use to initialize from.

    Returns
    -------
    instance: `~astropy.modeling.models`
        The initialized model.
        If there are any errors, the instance is returned
        uninitialized.
    """
    if x is None or y is None:
        return instance

    name = _get_model_name(instance)

    try:
        initializer = _initializers[name]()

        return initializer.initialize(instance, x, y)

    except KeyError:
        return instance
