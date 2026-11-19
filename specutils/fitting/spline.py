from __future__ import print_function, division, absolute_import

import numpy as np
from scipy import interpolate
import warnings

from astropy.modeling.core import FittableModel, Fittable1DModel, Model
from astropy.modeling.functional_models import Shift
from astropy.modeling.parameters import Parameter
from astropy.modeling.utils import poly_map_domain, comb
from astropy.modeling.fitting import _FitterMeta, fitter_unit_support
from astropy.utils import indent, check_broadcast
from astropy.units import Quantity

__all__ = []


class SplineModel(Fittable1DModel):
    """
    Wrapper around scipy.interpolate.splrep and splev
    
    Analogous to scipy.interpolate.UnivariateSpline() if knots unspecified,
    and scipy.interpolate.LSQUnivariateSpline if knots are specified
    
    There are two ways to make a spline model.
        1. you have the spline auto-determine knots from the data
        2. you specify the knots
    """

    def __init__(self, degree=3, smoothing=None, knots=None,
                 extrapolate_mode=0, *args, **kwargs):
        """
        Set up a spline model.

        degree: degree of the spline (default 3)
            In scipy fitpack, this is "k"

        smoothing (optional): smoothing value for automatically determining knots
            In scipy fitpack, this is "s"
            By default, uses s = len(w) (see scipy.interpolate.UnivariateSpline)

        knots (optional): spline knots (boundaries of piecewise polynomial)
            If not specified, will automatically determine knots based on
            degree + smoothing. The fit is identical to scipy.interpolate.UnivariateSpline.
            If specified, analogous to scipy.interpolate.LSQUnivariateSpline.

        extrapolate_mode (optional): how to deal with solution outside of interval.
            (see scipy.interpolate.splev)
            if 0 (default): return the extrapolated value
            if 1, return 0
            if 2, raise a ValueError
            if 3, return the boundary value
        """

        self._param_names = ()
        
        self._degree = degree
        self._smoothing = smoothing
        self._knots = self.verify_knots(knots)
        self.extrapolate_mode = extrapolate_mode

        # This is used to evaluate the spline. When None, raises an error when
        # trying to evaluate the spline.
        self._tck = None

        super().__init__(*args, **kwargs)
        
    def verify_knots(self, knots):
        """
        Basic knot array vetting. The goal of having this is to enable more
        useful error messages than scipy (if needed).
        """
        if knots is None:
            return None

        knots = np.array(knots)
        assert len(knots.shape) == 1, knots.shape
        knots = np.sort(knots)
        assert len(np.unique(knots)) == len(knots), knots

        return knots

    # Getters
    @property
    def degree(self):
        """ Spline degree (k in FITPACK) """
        return self._degree

    @property
    def smoothing(self):
        """ Spline smoothing (s in FITPACK) """
        return self._smoothing

    @property
    def knots(self):
        """ Spline knots (t in FITPACK) """
        return self._knots

    @property
    def coeffs(self):
        """ Spline coefficients (c in FITPACK) """
        if self._tck is not None:
            return self._tck[1]
        else:
            raise RuntimeError("SplineModel has not been fit yet.")

    # Setters
    def reset_model(self):
        """ Resets model so it needs to be refit to be valid """
        self._tck = None
        self._param_names = ()

    @degree.setter
    def degree(self, degree):
        """ Spline degree (k in FITPACK) """
        self._degree = degree
        self.reset_model()

    @smoothing.setter
    def smoothing(self, smoothing):
        """ Spline smoothing (s in FITPACK) """
        self._smoothing = smoothing
        self.reset_model()

    @knots.setter
    def knots(self, knots):
        """ Spline knots (t in FITPACK) """
        self._knots = self.verify_knots(knots)
        self.reset_model()

    def set_model_from_tck(self, tck):
        """
        Main way to update model
        Use output of scipy.interpolate.splrep
        """
        t, c, k = tck
        self.degree = k
        self.knots = t[k:-k]
        self._tck = tck
        self._param_names = self._generate_coeff_names()

    # Spline methods
    def derivative(self, n=1):
        if self._tck is None:
            raise RuntimeError("SplineModel has not been fit yet")
        else:
            ext = 1 if self.extrapolate_mode == 3 else self.extrapolate_mode
            new_tck = interpolate.fitpack.splder(self._tck, n)
            newmodel = SplineModel(degree=self.degree, smoothing=self.smoothing,
                                   knots=self.knots, extrapolate_mode=ext)
            newmodel.set_model_from_tck(new_tck)
            return newmodel

    def antiderivative(self, n=1):
        if self._tck is None:
            raise RuntimeError("SplineModel has not been fit yet.")
        else:
            new_tck = interpolate.fitpack.splantider(self._tck, n)
            newmodel = SplineModel(degree=self.degree, smoothing=self.smoothing,
                                   knots=self.knots, extrapolate_mode=self.extrapolate_mode)
            newmodel.set_model_from_tck(new_tck)
            return newmodel

    def integral(self, a, b):
        if self._tck is None:
            raise RuntimeError("SplineModel has not been fit yet.")
        else:
            t, c, k = self._tck
            return interpolate.dfitpack.splint(t, c, k, a, b)

    def derivatives(self, x):
        if self._tck is None:
            raise RuntimeError("SplineModel has not been fit yet.")
        else:
            t, c, k = self._tck
            d, ier = interpolate.dfitpack.spalde(t, c, k, x)
            if not ier == 0:
                raise ValueError("Error code returned by spalde: %s" % ier)
            return d

    def roots(self):
        if self._tck is None:
            raise RuntimeError("SplineModel has not been fit yet.")
        t, c, k = self._tck
        if k == 3:
            z, m, ier = interpolate.dfitpack.sproot(t, c)
            if not ier == 0:
                raise ValueError("Error code returned by spalde: %s" % ier)
            return z[:m]
        raise NotImplementedError('finding roots unsupported for '
                                  'non-cubic splines')

    def __call__(self, x, der=0):
        """
        Evaluate the model with the given inputs.
        der is passed to scipy.interpolate.splev
        """
        if self._tck is None:
            raise RuntimeError("SplineModel has not been fit yet.")

        return interpolate.splev(x, self._tck, der=der, ext=self.extrapolate_mode)

    # Stuff below here is stubs
    # TODO: fill out methods
    @property
    def param_names(self):
        """
        Coefficient names generated based on the model's knots and polynomial degree.
        Not Implemented
        """
        #raise NotImplementedError("SplineModel does not currently expose parameters")
        warnings.warn("SplineModel does not currently expose parameters\n"
                      "Will only work with SplineFitter")
        try:
            return self._param_names
        except AttributeError:
            return ()

    def _generate_coeff_names(self):
        names = []
        degree, Nknots = self._degree, len(self._knots)
        for i in range(Nknots):
            for j in range(degree+1):
                names.append("k{}_c{}".format(i,j))
        return tuple(names)

    def evaluate(self, *args, **kwargs):
        return self(*args, **kwargs)


class SplineFitter(metaclass=_FitterMeta):
    """
    Run a spline fit.
    """
    def __init__(self):
        self.fit_info = {"fp": None,
                         "ier": None,
                         "msg": None}
        super().__init__()

    def validate_model(self, model):
        if not isinstance(model, SplineModel):
            raise ValueError("model must be of type SplineModel (currently is {})".format(
                    type(model)))

    # TODO do something about units
    # @fitter_unit_support
    def __call__(self, model, x, y, w=None):
        """
        Fit a spline model to data.
        Internally uses scipy.interpolate.splrep.

        """

        self.validate_model(model)

        # Case (1): fit smoothing spline
        if model.knots is None:
            tck, fp, ier, msg = interpolate.splrep(x, y, w=w,
                                                   t=None,
                                                   k=model.degree,
                                                   s=model.smoothing,
                                                   task=0, full_output=True
                                                   )
        # Case (2): leastsq spline
        else:
            knots = model.knots
            ## TODO some sort of validation that the knots are internal, since
            ## this procedure automatically adds knots at the two endpoints
            tck, fp, ier, msg = interpolate.splrep(x, y, w=w,
                                                   t=knots,
                                                   k=model.degree,
                                                   s=model.smoothing,
                                                   task=-1, full_output=True
                                                   )

        model.set_model_from_tck(tck)
        self.fit_info.update({"fp": fp, "ier": ier, "msg": msg})
