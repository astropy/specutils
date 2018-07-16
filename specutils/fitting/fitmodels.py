from __future__ import division

import itertools
import operator

import numpy as np

from ..spectra.spectrum1d import Spectrum1D
from ..manipulation.utils import excise_regions
from astropy.modeling import fitting
import astropy.units as u


__all__ = ['fit_lines']


def fit_lines(spectrum, model, fitter=fitting.SLSQPLSQFitter(),
              exclude_regions=None, weights=None, window=None):
    """
    Fit the input models (initial conditions) to the spectrum.  Output will be
    the same models with the parameters set based on the fitting.

    spectrum, models -> models

    Parameters
    ----------
    spectrum : Spectrum1D
        The spectrum object overwhich the equivalent width will be calculated.

    models_initial_conds : ``astropy.modeling.models`` or list of ``astropy.modeling.models``
        The model or list of models that contain the initial guess.

    exclude_regions : list of 2-tuples
        List of regions to exclude in the fitting.

    weights : list  (NOT IMPLEMENTED YET)
        List of weights to define importance of fitting regions.

    Returns
    -------
    models : Compound model of ``astropy.modeling.models``
        A compound model of models with fitted parameters.

    Notes
    -----
       * Could add functionality to set the bounds in
         `model_initial_conds`` if they are not set.
       * The models in the list of `model_initial_conds` are added
          together and passed as a compound model to the
          ``astropy.modeling.fitter`` class instance.

    """

    #
    # If we are to exclude certain regions, then remove them.
    #

    if exclude_regions is not None:
        spectrum = excise_regions(spectrum, exclude_regions)

    #
    # Get the dispersion and flux to fit.
    #

    dispersion = spectrum.wavelength.value
    flux = spectrum.flux.value

    #
    # Make the model a list if not already
    #

    single_model_in = not isinstance(model, list)
    if not isinstance(model, list):
        model = [model]

    #
    # If a single model is passed in then just do that.
    #
    fitted_models = []

    for modeli, model_guess in enumerate(model):

        print('Window is {}'.format(window))

        #
        # Determine the window if it is not None.  There
        # are several options here:
        #   window = 4 * u.Angstrom -> Quantity
        #   window = (4*u.Angstrom, 6*u.Angstrom) -> tuple
        #   window = (4, 6)*u.Angstrom -> Quantity
        #

        #
        #  Determine the window if there is one
        #

        if window is not None and isinstance(window, list):
            model_window = window[modeli]
        elif window is not None:
            model_window = window
        else:
            model_window = None

        print('Model window is {}'.format(model_window))

        # In this case the window defines the area around the center of each model
        if model_window is not None and isinstance(model_window, (float, int)):
            center = model_guess.mean.value
            print('model_guaess is {}  model_window is {}'.format(center, model_window))
            indices = np.nonzero((dispersion >= center-model_window) & (dispersion < center+model_window))

            dispersion = dispersion[indices]
            flux = flux[indices]
            print('Going to fit to dispersion {}'.format(dispersion))

        # In this case the window is the start and end points of where we should fit
        elif model_window is not None and isinstance(model_window, tuple):
            indices = np.nonzero((dispersion >= model_window[0]) & (dispersion < model_window[1]))

            dispersion = dispersion[indices]
            flux = flux[indices]

        #
        # Do the fitting of spectrum to the model.
        #
        fit_model = fitter(model_guess, dispersion, flux)

        fitted_models.append(fit_model)

    #
    # If a single model was passed in as a paramter (and not a list)
    # then we are going to return just the model out.
    #
    if single_model_in:
        fitted_models = fitted_models[0]

    return fitted_models
