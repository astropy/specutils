from __future__ import division

from astropy import modeling
from astropy.modeling import models, fitting
from ..spectra import Spectrum1D

import numpy as np

import logging

__all__ = ['fit_continuum_generic']

def fit_continuum_generic(spectrum,
                          model=None, fitter=None,
                          sigma=3.0, sigma_lower=None, sigma_upper=None, iters=5,
                          exclude_regions=[],
                          full_output=False):
    """
    Fit a generic continuum model to a spectrum.
    
    The default algorithm is iterative sigma clipping
    
    Parameters
    ----------
    spectrum : `~specutils.Spectrum1D`
        The `~specutils.Spectrum1D` object to which a continuum model is fit

    model : `XXXX`
        The type of model to use for the continuum.
        astropy.modeling.models
        Must either be astropy.modeling.Fittable1DModel
        or the string "spline" (since this is not currently implemented)
        Default: models.Chebyshev1D(3)
        
    fitter : `XXXX`
        The type of fitter to use for the continuum.
        astropy.modeling.fitting
        Default: fitting.LevMarLSQFitter()
        
    sigma : float, optional
        The number of standard deviations to use for both lower and upper clipping limit.
        Defaults to 3.0
        
    sigma_lower : float or None, optional
        Number of standard deviations for lower bound clipping limit.
        If None (default), then `sigma` is used.
        
    sigma_upper : float or None, optional
        Number of standard deviations for upper bound clipping limit.
        If None (default), then `sigma` is used.
        
    iters : int or None, optional
        Number of iterations to perform sigma clipping.
        If None, clips until convergence achieved.
        Defaults to 5
        
    exclude_regions : list of tuples, optional
        A list of dispersion regions to exclude.
        Each tuple must be sorted.
        e.g. [(6555,6575)]

    full_output : bool, optional
        If True, return more information.
        Currently, just the model and the pixels-used boolean array
        
    Returns
    -------
    continuum_model : `XXXX`
        Output `XXXX` which is a model for the continuum

    Raises
    ------
    ValueError
       In the case that ``spectrum`` .... is not the correct type
    
    """
    
    ## Parameter checks
    if not isinstance(spectrum, Spectrum1D):
        raise ValueError('The spectrum parameter must be a Spectrum1D object')
    for exclude_region in exclude_regions:
        if len(exclude_region) != 2:
            raise ValueError('All exclusion regions must be of length 2')
        if exclude_region[0] >= exclude_region[1]:
            raise ValueError('All exclusion regions must be (low, high)')
    
    ## Set default model and fitter
    if model is None:
        logging.info("Using Chebyshev1D(3) as default continuum model")
        model = models.Chebyshev1D(3)
    if fitter is None:
        fitter = fitting.LevMarLSQFitter()
    if not isinstance(model, modeling.FittableModel):
        raise ValueError('The model parameter must be a astropy.modeling.FittableModel object')
    ## TODO this is waiting on a refactor in fitting to work
    #if not isinstance(fitter, fitting.Fitter):
    #    raise ValueError('The model parameter must be a astropy.modeling.fitting.Fitter object')

    ## Get input spectrum data
    x = spectrum.spectral_axis.value
    y = spectrum.flux.value
    
    ## Set up valid pixels mask
    ## Exclude non-finite values
    good = np.isfinite(y)
    ## Exclude regions
    for (excl1, excl2) in exclude_regions:
        good[np.logical_and(x > excl1, x < excl2)] = False
    
    ## Set up sigma clipping
    if sigma_lower is None: sigma_lower = sigma
    if sigma_upper is None: sigma_upper = sigma

    for i_iter in range(iters):
        logging.info("Iter {}: Fitting {}/{} pixels".format(i_iter, good.sum(), len(good)))
        ## Fit model
        ## TODO include data uncertainties
        new_model = fitter(model, x[good], y[good])
        
        ## Sigma clip
        difference = new_model(x) - y
        finite = np.isfinite(difference)
        sigma_difference = difference / np.std(difference[np.logical_and(good, finite)])
        good[sigma_difference > sigma_upper] = False
        good[sigma_difference < -sigma_lower] = False
        
        ## Update model iteratively: it is initialized at the previous fit's values
        #model = new_model

    model = new_model
    if full_output:
        return model, good
    return model
