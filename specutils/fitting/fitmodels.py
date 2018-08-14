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

    model_unitless, dispersion_unitless, flux_unitless = _strip_units_from_model(model, spectrum)

    #
    # Do the fitting of spectrum to the model.
    #

    print('model_unitless {}'.format(model_unitless))
    fit_model_unitless = fitter(model_unitless, dispersion_unitless, flux_unitless)

    #
    # Now add the units back onto the model....
    #

    fit_model = _add_units_to_model(fit_model_unitless, model)

    return fit_model


def _strip_units_from_model(model_in, spectrum):

    dispersion = spectrum.spectral_axis
    flux = spectrum.flux

    compound_model = model_in.n_submodels() > 1

    if not compound_model:
        model_in = [model_in]

    model_out = []
    for sub_model in model_in:

        new_sub_model = sub_model.__class__.copy(sub_model)

        for pn in new_sub_model.param_names:

            is_quantity = False
            a = getattr(sub_model, pn)
            if hasattr(a, 'quantity') and a.quantity is not None:
                is_quantity = True
                v = getattr(sub_model, pn).quantity.value * u.dimensionless_unscaled
            else:
                v = getattr(sub_model, pn).value

            setattr(new_sub_model, pn, v)
        model_out.append(new_sub_model)

    dispersion_out = dispersion.value
    if is_quantity:
        dispersion_out = dispersion_out * u.dimensionless_unscaled

    flux_out = flux.value
    if is_quantity:
        flux_out = flux_out * u.dimensionless_unscaled

    if compound_model:
        # TODO:  Wrong -- they may not be added together.
        model_out = functools.reduce(operator.add, model_out)
    else:
        model_out = model_out[0]

    return model_out, dispersion_out, flux_out


def _add_units_to_model(model_in, model_orig):

    compound_model = model_in.n_submodels() > 1
    if not compound_model:
        model_in = [model_in]
        model_orig = [model_orig]

    model_out = []
    for ii, m_in in enumerate(model_in):
        m_orig = model_orig[ii]

        new_sub_model = m_in.__class__.copy(m_in)

        for pn in new_sub_model.param_names:

            a = getattr(m_in, pn)
            if hasattr(a, 'quantity') and a.quantity is not None:
                v = getattr(m_in, pn).quantity.value * getattr(m_orig, pn).quantity.unit
            else:
                v = getattr(m_in, pn).value

            setattr(new_sub_model, pn, v)
        model_out.append(new_sub_model)

    if compound_model:
        # TODO:  Wrong -- they may not be added together.
        model_out = functools.reduce(operator.add, model_out)
    else:
        model_out = model_out[0]

    return model_out
