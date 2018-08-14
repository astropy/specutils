from __future__ import division

import functools
import operator

import numpy as np
import astropy.units as u

from ..spectra.spectrum1d import Spectrum1D
from ..manipulation.utils import excise_regions
from astropy.modeling import fitting


__all__ = ['fit_lines']


def fit_lines(spectrum, model, fitter=fitting.SimplexLSQFitter(),
              exclude_regions=None, weights=None, window=None):
    """
    Fit the input models (initial conditions) to the spectrum.  Output will be
    the same models with the parameters set based on the fitting.

    spectrum, models -> models

    Parameters
    ----------
    spectrum : Spectrum1D
        The spectrum object overwhich the equivalent width will be calculated.

    model: `~astropy.modeling.Model` or list of `~astropy.modeling.Model`
        The model or list of models that contain the initial guess.

    exclude_regions : list of `~specutils.SpectralRegion`
        List of regions to exclude in the fitting.

    weights : list  (NOT IMPLEMENTED YET)
        List of weights to define importance of fitting regions.

    Returns
    -------
    models : Compound model of `~astropy.modeling.Model`
        A compound model of models with fitted parameters.

    Notes
    -----
       * Could add functionality to set the bounds in
         ``model`` if they are not set.
       * The models in the list of ``model`` are added
          together and passed as a compound model to the
          `~astropy.modeling.fitting.Fitter` class instance.

    """

    #
    # If we are to exclude certain regions, then remove them.
    #

    if exclude_regions is not None:
        spectrum = excise_regions(spectrum, exclude_regions)

    #
    # Get the dispersion and flux to fit.
    #

    dispersion = spectrum.wavelength
    flux = spectrum.flux

    #
    # Make the model a list if not already
    #

    single_model_in = not isinstance(model, list)
    if single_model_in:
        model = [model]

    #
    # If a single model is passed in then just do that.
    #
    fitted_models = []

    for modeli, model_guess in enumerate(model):

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

        fit_model = _fit_lines(spectrum, model_guess, fitter,
                               exclude_regions, weights, model_window)

        fitted_models.append(fit_model)

    if single_model_in:
        fitted_models = fitted_models[0]

    return fitted_models


def _fit_lines(spectrum, model, fitter=fitting.SimplexLSQFitter(),
               exclude_regions=None, weights=None, window=None):
    """
    Fit the input model (initial conditions) to the spectrum.  Output will be
    the same model with the parameters set based on the fitting.

    spectrum, model -> model

    Parameters
    ----------
    spectrum : Spectrum1D
        The spectrum object overwhich the equivalent width will be calculated.

    model: `~astropy.modeling.Model`
        The model or that contain the initial guess.

    exclude_regions : list of `~specutils.SpectralRegion`
        List of regions to exclude in the fitting.

    weights : list  (NOT IMPLEMENTED YET)
        List of weights to define importance of fitting regions.

    Returns
    -------
    model : Compound model of `~astropy.modeling.Model`
        A compound model of models with fitted parameters.

    Notes
    -----
       * Could add functionality to set the bounds in ``model`` if they are not set.
       * The models in the list of ``model`` are added together and passed as a
          compound model to the `~astropy.modeling.fitting.Fitter` class instance.

    """

    #
    # If we are to exclude certain regions, then remove them.
    #

    if exclude_regions is not None:
        spectrum = excise_regions(spectrum, exclude_regions)

    dispersion = spectrum.spectral_axis
    flux = spectrum.flux

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

    # In this case the window defines the area around the center of each model
    if window is not None and isinstance(window, (float, int)):
        center = model.mean
        indices = np.nonzero((spectrum.spectral_axis >= center-window) & (spectrum.spectral_axis < center+window))

        dispersion = dispersion[indices]
        flux = flux[indices]

    # In this case the window is the start and end points of where we should fit
    elif window is not None and isinstance(window, tuple):
        indices = np.nonzero((dispersion >= window[0]) & (dispersion < window[1]))

        dispersion = dispersion[indices]
        flux = flux[indices]

    #
    # Compound models with units can not be fit.
    #
    # Convert the model initial guess to the spectral
    # units and then remove the units
    #

    model_unitless = _convert_and_remove_units(model, spectrum)

    #
    # Do the fitting of spectrum to the model.
    #

    fit_model_unitless = fitter(model_unitless, dispersion.value, flux.value)

    #
    # Now add the units back onto the model....
    #

    fit_model = _convert_and_add_units(fit_model_unitless, model, spectrum)

    return fit_model


def _convert_and_remove_units(model, spectrum):
    """
    This method converts the model's units to
    those of the spectrum and then outputs a new
    model with units stripped.

    # m.__class__(**{nm: getattr(m, nm).value for nm in m.param_names})
    """

    dispersion = spectrum.spectral_axis
    dispersion_unit = spectrum.spectral_axis.unit
    flux_unit = spectrum.flux.unit

    single_model_in = not hasattr(model, 'submodel_names')
    if single_model_in:
        model = [model]
        N_models = 1
    else:
        N_models = model.n_submodels()

    model_unitless = []
    for ii in range(N_models):
        m = model[ii]
        new_params = {}
        for param_name in m.param_names:
            quantity = getattr(m, param_name).quantity

            if quantity is not None:

                if quantity.unit.is_equivalent(dispersion_unit, equivalencies=u.equivalencies.spectral()):
                    quantity = quantity.to(dispersion_unit, equivalencies=u.equivalencies.spectral())

                elif quantity.unit.is_equivalent(flux_unit, equivalencies=u.equivalencies.spectral_density(dispersion)):
                    quantity = quantity.to(flux_unit, equivalencies=u.equivalencies.spectral_density(dispersion))

                else:
                    raise Exception('Conversion from unit type {}'.format(quantity.unit))

                new_params[param_name] = quantity.value
            else:
                new_params[param_name] = getattr(m, param_name).value

        # Now that all the parameters have been cleaned up
        # create the new model class
        model_guess = m.deepcopy()
        [setattr(model_guess, pn, new_params[pn]) for pn in model_guess.param_names]

        model_unitless.append(model_guess)

    if single_model_in:
        return model_unitless[0]
    else:
        return functools.reduce(operator.add, model_unitless)


def _convert_and_add_units(model, model_init, spectrum):
    """
    This method converts the model's units to
    those of the spectrum and then outputs a new
    model with units stripped.

    # m.__class__(**{nm: getattr(m, nm).value for nm in m.param_names})
    """

    single_model_in = not hasattr(model, 'submodel_names')
    if single_model_in:
        model = [model]
        model_init = [model_init]
        N_models = 1
    else:
        N_models = model.n_submodels()

    model_unitless = []

    for ii in range(N_models):
        m = model[ii]
        m_init = model_init[ii]

        new_params = {}
        for param_name in m.param_names:
            quantity = getattr(m, param_name).value

            if getattr(m_init, param_name).quantity is not None:
                new_params[param_name] = quantity * getattr(m_init, param_name).quantity.unit
            else:
                new_params[param_name] = quantity

        # Now that all the parameters have been cleaned up
        # create the new model class
        #model_out = m.__class__(**new_params)
        model_out = m.deepcopy()
        [setattr(model_out, pn, u.Quantity(new_params[pn])) for pn in model_out.param_names]

        model_unitless.append(model_out)

    if single_model_in:
        return model_unitless[0]
    else:
        return functools.reduce(operator.add, model_unitless)
