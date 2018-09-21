from __future__ import division

import operator

import numpy as np
import astropy.units as u

from ..manipulation.utils import excise_regions
from ..utils import QuantityModel
from astropy.modeling import fitting, Model, models


__all__ = ['fit_lines']


def fit_lines(spectrum, model, fitter=fitting.LevMarLSQFitter(),
              exclude_regions=None, weights=None, window=None):
    """
    Fit the input models to the spectrum. The parameter values of the
    input models will be used as the initial conditions for the fit.

    Parameters
    ----------
    spectrum : Spectrum1D
        The spectrum object over which the equivalent width will be calculated.

    model: `~astropy.modeling.Model` or list of `~astropy.modeling.Model`
        The model or list of models that contain the initial guess.

    exclude_regions : list of `~specutils.SpectralRegion`
        List of regions to exclude in the fitting.

    weights : list  (NOT IMPLEMENTED YET)
        List of weights to define importance of fitting regions.

    window : `~specutils.SpectralRegion` or list of `~specutils.SpectralRegion`
        Regions of the spectrum to use in the fitting. If None, then the
        whole spectrum will be used in the fitting.

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


def _fit_lines(spectrum, model, fitter=fitting.LevMarLSQFitter(),
               exclude_regions=None, weights=None, window=None):
    """
    Fit the input model (initial conditions) to the spectrum.  Output will be
    the same model with the parameters set based on the fitting.

    spectrum, model -> model

    Parameters
    ----------
    spectrum : Spectrum1D
        The spectrum object over which the equivalent width will be calculated.

    model: `~astropy.modeling.Model`
        The model or that contain the initial guess.

    exclude_regions : list of `~specutils.SpectralRegion`
        List of regions to exclude in the fitting.

    weights : list  (NOT IMPLEMENTED YET)
        List of weights to define importance of fitting regions.

    window : `~specutils.SpectralRegion` or list of `~specutils.SpectralRegion`
        Regions of the spectrum to use in the fitting. If None, then the
        whole spectrum will be used in the fitting.

    Returns
    -------
    model : Compound model of `~astropy.modeling.Model`
        A compound model of models with fitted parameters.

    Notes
    -----
       * Could add functionality to set the bounds in ``model`` if they are not set.

    """

    if weights is not None:
        raise NotImplementedError('Weights are not yet implemented.')

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

    fit_model_unitless = fitter(model_unitless, dispersion_unitless, flux_unitless)

    #
    # Now add the units back onto the model....
    #

    fit_model = _add_units_to_model(fit_model_unitless, model, spectrum)

    return fit_model


def _strip_units_from_model(model_in, spectrum):
    """
    This method strips the units from the model, so the result can
    be passed to the fitting routine. This is necessary as CoumpoundModel
    with units does not work in the fitters.

    Note:  When CompoundModel with units works in the fitters this method
           can be removed.

    Note: This assumes there are two types of models, those that are
          based on `~astropy.modeling.models.PolynomialModel` and therefore
          require the ``degree`` parameter when instantiating the class, and
          "everything else" that does not require an "extra" parameter for
          class instantiation.
    """

    #
    # Get the dispersion and flux information from the spectrum
    #

    dispersion = spectrum.spectral_axis
    dispersion_unit = spectrum.spectral_axis.unit

    flux = spectrum.flux
    flux_unit = spectrum.flux.unit

    #
    # Determine if a compound model
    #

    compound_model = model_in.n_submodels() > 1

    if not compound_model:
        # For this we are going to just make it a list so that we
        # can use the looping structure below.
        model_in = [model_in]
    else:
        # If it is a compound model then we are going to create the RPN
        # representation of it which is a list that contains either astropy
        # models or string representations of operators (e.g., '+' or '*').
        model_in = [c.value for c in model_in._tree.traverse_postorder()]

    #
    # Run through each model in the list or compound model
    #

    model_out_stack = []
    for sub_model in model_in:

        #
        # If it is an operator put onto the stack and move on...
        #

        if not isinstance(sub_model, Model):
            model_out_stack.append(sub_model)
            continue

        #
        # Make a new instance of the class.
        #

        if isinstance(sub_model, models.PolynomialModel):
            new_sub_model = sub_model.__class__(sub_model.degree)
        else:
            new_sub_model = sub_model.__class__()

        #
        # Now for each parameter in the model determine if a dispersion or
        # flux type of unit, then convert to spectrum units and then get the value.
        #

        for pn in new_sub_model.param_names:

            # This could be a Quantity or Parameter
            a = getattr(sub_model, pn)

            if hasattr(a, 'quantity') and a.quantity is not None:
                q = getattr(sub_model, pn).quantity

                #
                # Convert the quantity to the spectrum's units, and then we will use
                # the *value* of it in the new unitless-model.
                #

                if q.unit.is_equivalent(dispersion_unit, equivalencies=u.equivalencies.spectral()):
                    quantity = q.to(dispersion_unit, equivalencies=u.equivalencies.spectral())

                elif q.unit.is_equivalent(flux_unit, equivalencies=u.equivalencies.spectral_density(dispersion)):
                    quantity = q.to(flux_unit, equivalencies=u.equivalencies.spectral_density(dispersion))

                # The value must be a quantity in order to use setattr (below)
                # as the setter requires a Quantity on the rhs if the parameter
                # is already a value.
                v = quantity.value
            else:
                v = getattr(sub_model, pn).value

            #
            # Add this information for the parameter name into the
            # new sub model.
            #

            setattr(new_sub_model, pn, v)

        # The new model now has unitless information in it but has
        # been converted to spectral unit scale.
        model_out_stack.append(new_sub_model)

    # If a compound model we need to re-create it, otherwise
    # it is a single model and we just get the first one (as
    # there is only one).
    if compound_model:
        model_out = _combine_postfix(model_out_stack)
    else:
        model_out = model_out_stack[0]

    return model_out, dispersion.value, flux.value


def _add_units_to_model(model_in, model_orig, spectrum):
    """
    This method adds the units to the model based on the units of the
    model passed in.  This is necessary as CoumpoundModel
    with units does not work in the fitters.

    Note:  When CompoundModel with units works in the fitters this method
           can be removed.

    Note: This assumes there are two types of models, those that are
          based on `~astropy.modeling.models.PolynomialModel` and therefore
          require the ``degree`` parameter when instantiating the class, and
          "everything else" that does not require an "extra" parameter for
          class instantiation.
    """

    dispersion = spectrum.spectral_axis

    #
    # If not a compound model, then make a single element
    # list so we can use the for loop below.
    #

    compound_model = model_in.n_submodels() > 1
    if not compound_model:
        model_in_list = [model_in]
        model_orig_list = [model_orig]
    else:
        compound_model_in = model_in

        model_in_list = [c.value for c in model_in._tree.traverse_postorder()]
        model_orig_list = [c.value for c in model_orig._tree.traverse_postorder()]

    model_out_stack = []
    model_index = 0

    #
    # For each model in the list we will convert the values back to
    # the original (sub-)model units.
    #

    for ii, m_in in enumerate(model_in_list):

        #
        # If an operator (ie not Model) then we'll just add
        # to the stack and evaluate at the end.
        #

        if not isinstance(m_in, Model):
            model_out_stack.append(m_in)
            continue

        #
        # Get the corresponding *original* sub-model that
        # will match the current sub-model. From this we will
        # grab the units to apply.
        #

        m_orig = model_orig_list[ii]

        #
        # Make the new sub-model.
        #

        if isinstance(m_in, models.PolynomialModel):
            new_sub_model = m_in.__class__(m_in.degree)
        else:
            new_sub_model = m_in.__class__()

        #
        # Convert the model values from the spectrum units back to the
        # original model units.
        #

        for pi, pn in enumerate(new_sub_model.param_names):

            #
            # Get the parameter from the original model and unit-less model.
            #

            m_orig_param = getattr(m_orig, pn)
            m_in_param = getattr(m_in, pn)

            if hasattr(m_orig_param, 'quantity') and m_orig_param.quantity is not None:

                m_orig_param_quantity = m_orig_param.quantity

                #
                # If a spectral dispersion type of unit...
                #
                if m_orig_param_quantity.unit.is_equivalent(spectrum.spectral_axis.unit,
                                                            equivalencies=u.equivalencies.spectral()):

                    # If it is a compound model, then we need to get the value from the
                    # actual compound model as the tree is not updated in the fitting
                    if compound_model:
                        current_value = getattr(compound_model_in, '{}_{}'.format(pn, model_index)).value *\
                                        spectrum.spectral_axis.unit
                    else:
                        current_value = m_in_param.value * spectrum.spectral_axis.unit

                    v = current_value.to(m_orig_param_quantity.unit, equivalencies=u.equivalencies.spectral())

                #
                # If a spectral density type of unit...
                #
                elif m_orig_param_quantity.unit.is_equivalent(spectrum.flux.unit,
                                                              equivalencies=u.equivalencies.spectral_density(dispersion)):
                    # If it is a compound model, then we need to get the value from the
                    # actual compound model as the tree is not updated in the fitting
                    if compound_model:
                        current_value = getattr(compound_model_in, '{}_{}'.format(pn, model_index)).value *\
                                        spectrum.flux.unit
                    else:
                        current_value = m_in_param.value * spectrum.flux.unit

                    v = current_value.to(m_orig_param_quantity.unit,
                                         equivalencies=u.equivalencies.spectral_density(dispersion))

            else:
                v = getattr(m_in, pn).value

            #
            # Set the parameter value into the new sub-model.
            #

            setattr(new_sub_model, pn, v)

        #
        # Add the new unit-filled model onto the stack.
        #

        model_out_stack.append(new_sub_model)

        model_index += 1

    #
    # Create the output model which is either the evaulation
    # of the RPN representation of the model (if a compound model)
    # or just the first element if a non-compound model.
    #

    if compound_model:
        model_out = _combine_postfix(model_out_stack)
    else:
        model_out = model_out_stack[0]

    # If the first parameter is not a Quantity, then at this point we will assume
    # none of them are. (It would be inconsistent for fitting to have a model that
    # has some parameters as Quantities and some values).
    if getattr(model_orig, model_orig.param_names[0]).unit is None:
        model_out = QuantityModel(model_out, spectrum.spectral_axis.unit, spectrum.flux.unit)

    return model_out


def _combine_postfix(equation):
    """
    Given a Python list in post order (RPN) of an equation, convert/apply the operations to evaluate.
    The list order is the same as what is output from ``model._tree.traverse_postorder()``.

    Structure modified from https://codereview.stackexchange.com/questions/79795/reverse-polish-notation-calculator-in-python
    """

    ops = {'+': operator.add,
           '-': operator.sub,
           '*': operator.mul,
           '/': operator.truediv,
           '^': operator.pow,
           '**': operator.pow}

    stack = []
    result = 0
    for i in equation:
        if isinstance(i, Model):
            stack.insert(0, i)
        else:
            if len(stack) < 2:
                print('Error: insufficient values in expression')
                break
            else:
                n1 = stack.pop(1)
                n2 = stack.pop(0)
                result = ops[i](n1, n2)
                stack.insert(0, result)
    return result
