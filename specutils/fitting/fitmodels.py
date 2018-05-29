from __future__ import division

import itertools
import operator

import numpy as np
from astropy.modeling import models, fitting


__all__ = ['fit_models', 'fit_models_simple']


def fit_models(spectrum, model_initial_conds, fit_models_type='simple',
               window=None, weights=None, *args, **kwargs):
    """
    Does a naive equivalent width measures on the spectrum object.

    Parameters
    ----------
    spectrum : Spectrum1D
        The spectrum object overwhich the equivalent width will be calculated.
    model_initial_conds : list of ``astropy.modeling.models``
        The list of models that contain the initial guess.
    fit_models_type: str
        String representation of fit method to use as defined by the dict fit_models_types.
    window : tuple of wavelengths 
        Start and end wavelengths used for fitting.
    weights : list  (NOT IMPLEMENTED YET)
        List of weights to define importance of fitting regions.

    Returns
    -------
    models : list of ``astropy.modeling.models``
        The list of models that contain the fitted model.

    """

    #
    #  Define the fit line methods that are available.
    #     str description  ->  fit method
    #

    fit_models_types = {
        'simple': fit_models_simple
    }

    if fit_models_type in fit_models_types.keys():
        return fit_models_types[fit_models_type](spectrum, model_initial_conds, 
                                                 *args, **kwargs)
    else:
        raise ValueError('Fit line type {} is not one of implmented {}'.format(
            fit_models_type, fit_models_types.keys()))


def fit_models_simple(spectrum, model_initial_conds, fitter=fitting.SLSQPLSQFitter, window=None, weights=None):
    """
    Does a naive equivalent width measures on the spectrum object.

    Parameters
    ----------
    spectrum : Spectrum1D
        The spectrum object overwhich the equivalent width will be calculated.
    models : list of ``astropy.modeling.models``
        The list of models that contain the initial guess.
    window : tuple of wavelengths 
        Start and end wavelengths used for fitting.
    weights : list  (NOT IMPLEMENTED YET)
        List of weights to define importance of fitting regions.

    Returns
    -------
    models : list of ``astropy.modeling.models``
        The list of models that contain the fitted model.

    """

    print(fitter)

    #
    # Now to fit the data create a new superposition with initial
    # guesses for the parameters:
    #

    *_, compound_model = itertools.accumulate(model_initial_conds, operator.add)

    #
    # Do the fitting of spectrum to the model.
    #

    fitter_func = fitter()
    fitted_models = fitter_func(compound_model, spectrum.wavelength.value, spectrum.flux.value)

    #
    # Split up the combined model in the component form.
    #

    return fitted_models
