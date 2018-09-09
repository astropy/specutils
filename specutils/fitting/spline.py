from __future__ import print_function, division, absolute_import

import numpy as np
from scipy import interpolate

from astropy.modeling.core import FittableModel, Model
from astropy.modeling.functional_models import Shift
from astropy.modeling.parameters import Parameter
from astropy.modeling.utils import poly_map_domain, comb
from astropy.modeling.fitting import _FitterMeta, fitter_unit_support
from astropy.utils import indent, check_broadcast
from astropy.units import Quantity


__all__ = []

class SplineModel(FittableModel):
    """
    Wrapper around scipy.interpolate.splrep and splev
    
    Analogous to scipy.interpolate.UnivariateSpline() if knots unspecified,
    and scipy.interpolate.LSQUnivariateSpline if knots are specified
    
    There are two ways to make a spline model.
    (1) you have the spline auto-determine knots from the data
    (2) you specify the knots
    
    """
    
    linear = False # I think? I have no idea?
    col_fit_deriv = False # Not sure what this is
    
    def __init__(self, degree=3, smoothing=None, knots=None, extrapolate_mode=0):
        """
        Set up a spline model.
        
        degree: degree of the spline (default 3)
            In scipy fitpack, this is "k"
        
        smoothing (optional): smoothing value for automatically determining knots
            In scipy fitpack, this is "s"
            By default, uses a 
        
        knots (optional): spline knots (boundaries of piecewise polynomial)
            If not specified, will automatically determine knots based on
            degree + smoothing
            
        extrapolate_mode (optional): how to deal with solution outside of interval.
            (see scipy.interpolate.splev)
            if 0 (default): return the extrapolated value
            if 1, return 0
            if 2, raise a ValueError
            if 3, return the boundary value
        """
        self._degree = degree
        self._smoothing = smoothing
        self._knots = self.verify_knots(knots)
        self.extrapolate_mode = extrapolate_mode
        
        ## This is used to evaluate the spline
        ## When None, raises an error when trying to evaluate the spline
        self._tck = None
        
        self._param_names = ()
        
    def verify_knots(self, knots):
        """
        Basic knot array vetting.
        The goal of having this is to enable more useful error messages
        than scipy (if needed).
        """
        if knots is None: return None
        knots = np.array(knots)
        assert len(knots.shape) == 1, knots.shape
        knots = np.sort(knots)
        assert len(np.unique(knots)) == len(knots), knots
        return knots
    
    ############
    ## Getters
    ############
    def get_degree(self):
        """ Spline degree (k in FITPACK) """
        return self._degree
    def get_smoothing(self):
        """ Spline smoothing (s in FITPACK) """
        return self._smoothing
    def get_knots(self):
        """ Spline knots (t in FITPACK) """
        return self._knots
    def get_coeffs(self):
        """ Spline coefficients (c in FITPACK) """
        if self._tck is not None:
            return self._tck[1]
        else:
            raise RuntimeError("SplineModel has not been fit yet")
    
    ############
    ## Spline methods: not tested at all
    ############
    def derivative(self, n=1):
        if self._tck is None:
            raise RuntimeError("SplineModel has not been fit yet")
        else:
            t, c, k = self._tck
            return scipy.interpolate.BSpline.construct_fast(
                t,c,k,extrapolate=(self.extrapolate_mode==0)).derivative(n)
    def antiderivative(self, n=1):
        if self._tck is None:
            raise RuntimeError("SplineModel has not been fit yet")
        else:
            t, c, k = self._tck
            return scipy.interpolate.BSpline.construct_fast(
                t,c,k,extrapolate=(self.extrapolate_mode==0)).antiderivative(n)
    def integral(self, a, b):
        if self._tck is None:
            raise RuntimeError("SplineModel has not been fit yet")
        else:
            t, c, k = self._tck
            return scipy.interpolate.BSpline.construct_fast(
                t,c,k,extrapolate=(self.extrapolate_mode==0)).integral(a,b)
    def derivatives(self, x):
        raise NotImplementedError
    def roots(self):
        raise NotImplementedError
    
    ############
    ## Setters: not really implemented or tested
    ############
    def reset_model(self):
        """ Resets model so it needs to be refit to be valid """
        self._tck = None
    def set_degree(self, degree):
        """ Spline degree (k in FITPACK) """
        raise NotImplementedError
        self._degree = degree
        self.reset_model()
    def set_smoothing(self, smoothing):
        """ Spline smoothing (s in FITPACK) """
        raise NotImplementedError
        self._smoothing = smoothing
        self.reset_model()
    def set_knots(self, knots):
        """ Spline knots (t in FITPACK) """
        raise NotImplementedError
        self._knots = self.verify_knots(knots)
        self.reset_model()
    
    def set_model_from_tck(self, tck):
        """
        Use output of scipy.interpolate.splrep
        """
        self._tck = tck

    def __call__(self, x, der=0):
        """
        Evaluate the model with the given inputs.
        der is passed to scipy.interpolate.splev
        """
        if self._tck is None:
            raise RuntimeError("SplineModel has not been fit yet")
        return interpolate.splev(x, self._tck, der=der, ext=self.extrapolate_mode)
    
    ####################################
    ######### Stuff below here is stubs
    @property
    def param_names(self):
        """
        Coefficient names generated based on the model's knots and polynomial degree.
        Not Implemented
        """
        raise NotImplementedError("SplineModel does not currently expose parameters")
        return self._param_names

    #def __getattr__(self, attr):
    #    """
    #    Fails right now. Future code:
    #    # From astropy.modeling.polynomial.PolynomialBase
    #    if self._param_names and attr in self._param_names:
    #        return Parameter(attr, default=0.0, model=self)
    #    raise AttributeError(attr)
    #    """
    #    raise NotImplementedError("SplineModel does not currently expose parameters")

    #def __setattr__(self, attr, value):
    #    """
    #    Fails right now. Future code:
    #    # From astropy.modeling.polynomial.PolynomialBase
    #    if attr[0] != '_' and self._param_names and attr in self._param_names:
    #        param = Parameter(attr, default=0.0, model=self)
    #        param.__set__(self, value)
    #    else:
    #        super().__setattr__(attr, value)
    #    """
    #    raise NotImplementedError("SplineModel does not currently expose parameters")

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
    
    ## TODO do something about units
    #@fitter_unit_support
    def __call__(self, model, x, y, w=None):
        """
        Fit a spline model to data.
        Internally uses scipy.interpolate.splrep.
        
        """
        
        self.validate_model(model)
        
        ## Case (1): fit smoothing spline
        if model.get_knots() is None:
            tck, fp, ier, msg = interpolate.splrep(x, y, w=w,
                                                   t=None,
                                                   k=model.get_degree(), 
                                                   s=model.get_smoothing(),
                                                   task=0, full_output=True
                                                   )
        ## Case (2): leastsq spline
        else:
            knots = model.get_knots()
            ## TODO some sort of validation that the knots are internal, since
            ## this procedure automatically adds knots at the two endpoints
            tck, fp, ier, msg = interpolate.splrep(x, y, w=w,
                                                   t=knots,
                                                   k=model.get_degree(), 
                                                   s=model.get_smoothing(),
                                                   task=-1, full_output=True
                                                   )
        
        model.set_model_from_tck(tck)
        self.fit_info.update({"fp":fp, "ier":ier, "msg":msg})
    
