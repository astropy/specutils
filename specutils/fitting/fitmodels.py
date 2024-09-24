import itertools
import operator
import warnings

import numpy as np
from astropy import units as u
from astropy.nddata import StdDevUncertainty
from astropy.modeling import fitting, Model, models
from astropy.table import QTable
from scipy.signal import convolve

from ..spectra.spectral_region import SpectralRegion
from ..spectra.spectrum1d import Spectrum1D
from ..utils import QuantityModel
from ..analysis import fwhm, gaussian_sigma_width, centroid, warn_continuum_below_threshold
from ..manipulation import extract_region
from ..manipulation.utils import excise_regions

__all__ = ['find_lines_threshold', 'find_lines_derivative', 'fit_lines',
           'estimate_line_parameters']

# Define the initial estimators. This are the default methods to use to
# estimate astropy model parameters. This is based on only a small subset of
# the astropy models but it was determined that this is a decent start as most
# fitting will probably use one of these.
#
# Each method list must take a Spectrum1D object and should return a Quantity.
_parameter_estimators = {
    'Gaussian1D': {
        'amplitude': lambda s: max(s.flux),
        'mean': lambda s, region: centroid(s, regions=region),
        'stddev': lambda s, region: gaussian_sigma_width(s, regions=region)
    },
    'Lorentz1D': {
        'amplitude': lambda s: max(s.flux),
        'x_0': lambda s, region: centroid(s, regions=region),
        'fwhm': lambda s, region: fwhm(s, regions=region)
    },
    'Voigt1D': {
        'x_0': lambda s, region: centroid(s, regions=region),
        'amplitude_L': lambda s: max(s.flux),
        'fwhm_L': lambda s, region: fwhm(s, regions=region) / np.sqrt(2),
        'fwhm_G': lambda s, region: fwhm(s, regions=region) / np.sqrt(2)
    }
}


def _set_parameter_estimators(model):
    """
    Helper method used in method below.
    """
    if model.__class__.__name__ in _parameter_estimators:
        model_pars = _parameter_estimators[model.__class__.__name__]
        for name in model.param_names:
            par = getattr(model, name)
            setattr(par, "estimator", model_pars[name])

    return model


def estimate_line_parameters(spectrum, model, region=None):
    """
    The input ``model`` parameters will be estimated from the input
    ``spectrum``. The ``model`` can be specified with default parameters, for
    example ``Gaussian1D()``.

    Parameters
    ----------
    spectrum : `~specutils.Spectrum1D`
        The spectrum object from which we will estimate the model parameters.

    model : `~astropy.modeling.Model`
        Model for which we want to estimate parameters from the spectrum.

    Returns
    -------
    model : `~astropy.modeling.Model`
        Model with parameters estimated.
    """

    model = _set_parameter_estimators(model)
    # Estimate the parameters based on the estimators already
    # attached to the model
    for name in model.param_names:
        par = getattr(model, name)
        try:
            estimator = getattr(par, "estimator")
            if name[:9] == "amplitude":
                setattr(model, name, estimator(spectrum))
            else:
                setattr(model, name, estimator(spectrum, region))
        except AttributeError:
            raise Exception('No method to estimate parameter {}'.format(name))

    return model


def _consecutive(data, stepsize=1):
    return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)


@warn_continuum_below_threshold(threshold=0.01)
def find_lines_threshold(spectrum, noise_factor=1):
    """
    Find the emission and absorption lines in a spectrum. The method
    here is based on deviations larger than the spectrum's uncertainty
    (converted to standard deviation if in another representation) by the
    ``noise_factor``.

    This method only works with continuum-subtracted spectra and the
    uncertainty must be defined on the spectrum. To add the uncertainty,
    one could use `~specutils.manipulation.noise_region_uncertainty` to add
    the uncertainty.

    Parameters
    ----------
    spectrum : `~specutils.Spectrum1D`
        The spectrum object in which the lines will be found.

    noise_factor : float
       ``noise_factor`` multiplied by the spectrum's``uncertainty``, used for
        thresholding.

    Returns
    -------
    qtable: `~astropy.table.QTable`
        Table of emission and absorption lines. Line center (``line_center``),
        line type (``line_type``) and index of line center
        (``line_center_index``) are stored for each line.
    """

    # Threshold based on noise estimate and factor.
    uncertainty = spectrum.uncertainty
    uncert_val = 0 if uncertainty is None else uncertainty.represent_as(StdDevUncertainty).array

    inds = np.where(np.abs(spectrum.flux) > (noise_factor * uncert_val) *
                    spectrum.flux.unit)[0]
    pos_inds = inds[spectrum.flux.value[inds] > 0]
    line_inds_grouped = _consecutive(pos_inds, stepsize=1)

    if len(line_inds_grouped[0]) > 0:
        emission_inds = [inds[np.argmax(spectrum.flux.value[inds])]
                         for inds in line_inds_grouped]
    else:
        emission_inds = []

    #
    # Find the absorption lines
    #

    neg_inds = inds[spectrum.flux.value[inds] < 0]
    line_inds_grouped = _consecutive(neg_inds, stepsize=1)

    if len(line_inds_grouped[0]) > 0:
        absorption_inds = [inds[np.argmin(spectrum.flux.value[inds])]
                           for inds in line_inds_grouped]
    else:
        absorption_inds = []

    return _generate_line_list_table(spectrum, emission_inds, absorption_inds)


@warn_continuum_below_threshold(threshold=0.01)
def find_lines_derivative(spectrum, flux_threshold=None):
    """
    Find the emission and absorption lines in a spectrum. The method
    here is based on finding the zero crossings in the derivative
    of the spectrum.

    Parameters
    ----------
    spectrum : `~specutils.Spectrum1D`
        The spectrum object over which the equivalent width will be calculated.
    flux_threshold : float, `~astropy.units.Quantity` or None
        The threshold a pixel must be above to be considered part of a line. If
        a float, will assume the same units as ``spectrum.flux``. This
        threshold is above and beyond the derivative searching step. Default
        is None so no thresholding. The threshold is positive for emission
        lines and negative for absorption lines.

    Returns
    -------
    qtable: `~astropy.table.QTable`
        Table of emission and absorption lines. Line center (``line_center``),
        line type (``line_type``) and index of line center
        (``line_center_index``) are stored for each line.
    """

    # Take the derivative to find the zero crossings which correspond to
    # the peaks (positive or negative)
    kernel = [1, 0, -1]
    dY = convolve(spectrum.flux, kernel, 'valid')

    # Use sign flipping to determine direction of change
    S = np.sign(dY)
    ddS = convolve(S, kernel, 'valid')

    # Add units if needed.
    if flux_threshold is not None and isinstance(flux_threshold, (int, float)):
        flux_threshold = float(flux_threshold) * spectrum.flux.unit

    #
    # Emmision lines
    #

    # Find all the indices that appear to be part of a +ve peak
    candidates = np.where(dY > 0)[0] + (len(kernel) - 1)
    line_inds = sorted(set(candidates).intersection(np.where(ddS == -2)[0] + 1))

    if flux_threshold is not None:
        line_inds = np.array(line_inds)[spectrum.flux[line_inds] > flux_threshold]

    # Now group them and find the max highest point.
    line_inds_grouped = _consecutive(line_inds, stepsize=1)

    if len(line_inds_grouped[0]) > 0:
        emission_inds = [inds[np.argmax(spectrum.flux[inds])]
                         for inds in line_inds_grouped]
    else:
        emission_inds = []

    #
    # Absorption lines
    #

    # Find all the indices that appear to be part of a -ve peak
    candidates = np.where(dY < 0)[0] + (len(kernel) - 1)
    line_inds = sorted(set(candidates).intersection(np.where(ddS == 2)[0] + 1))

    if flux_threshold is not None:
        line_inds = np.array(line_inds)[spectrum.flux[line_inds] < -flux_threshold]

    # Now group them and find the max highest point.
    line_inds_grouped = _consecutive(line_inds, stepsize=1)

    if len(line_inds_grouped[0]) > 0:
        absorption_inds = [inds[np.argmin(spectrum.flux[inds])] for inds in
                           line_inds_grouped]
    else:
        absorption_inds = []

    return _generate_line_list_table(spectrum, emission_inds, absorption_inds)


def _generate_line_list_table(spectrum, emission_inds, absorption_inds):
    qtable = QTable()
    qtable['line_center'] = list(
        itertools.chain(
            *[spectrum.spectral_axis.value[emission_inds],
              spectrum.spectral_axis.value[absorption_inds]]
        )) * spectrum.spectral_axis.unit
    qtable['line_type'] = ['emission'] * len(emission_inds) + \
                          ['absorption'] * len(absorption_inds)
    qtable['line_center_index'] = list(
        itertools.chain(
            *[emission_inds, absorption_inds]))

    return qtable


def fit_lines(spectrum, model, fitter=fitting.LevMarLSQFitter(calc_uncertainties=True),
              exclude_regions=None, exclude_region_upper_bounds=False, weights=None,
              window=None, get_fit_info=False, **kwargs):
    """
    Fit the input models to the spectrum. The parameter values of the
    input models will be used as the initial conditions for the fit.

    Parameters
    ----------
    spectrum : `~specutils.Spectrum1D`
        The spectrum object over which the equivalent width will be calculated.
    model: `~astropy.modeling.Model` or list of `~astropy.modeling.Model`
        The model or list of models that contain the initial guess.
    fitter : `~astropy.modeling.fitting.Fitter`, optional
        Fitter instance to be used when fitting model to spectrum.
    exclude_regions : list of `~specutils.SpectralRegion`
        List of regions to exclude in the fitting.
    exclude_region_upper_bounds : bool, optional
        Set to True to override the default pythonic slicing on the excluded regions (inclusive
        lower bound, exclusive upper bound) and remove spectral values less than or equal to
        the upper bound of the region rather than strictly less than the upper bound.
    weights : array-like or 'unc', optional
        If 'unc', the uncertainties from the spectrum object are used to
        to calculate the weights. If array-like, represents the weights to
        use in the fitting.  Note that if a mask is present on the spectrum, it
        will be applied to the ``weights`` as it would be to the spectrum
        itself.
    window : `~specutils.SpectralRegion`, `~astropy.units.Quantity`, or list of either
        Regions of the spectrum to use in the fitting. If None, then the
        whole spectrum will be used in the fitting. If a single `~astropy.units.Quantity`
        is input, it will be used as the width of the region around the model mean.
    get_fit_info : bool, optional
        Flag to return the ``fit_info`` from the underlying scipy optimizer used
        in the fitting. If True, the returned model will have a ``fit_info``
        key populated in its ``meta`` dictionary.
    Additional keyword arguments are passed directly into the call to the
    ``fitter``.

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
        spectrum = excise_regions(spectrum, exclude_regions,
                                  inclusive_upper=exclude_region_upper_bounds)

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

        #
        # Check to see if the model has units. If it does not
        # have units then we are going to ignore them.
        #

        ignore_units = getattr(model_guess, model_guess.param_names[0]).unit is None

        fit_model = _fit_lines(
            spectrum,
            model_guess,
            fitter=fitter,
            exclude_regions=exclude_regions,
            exclude_region_upper_bounds=exclude_region_upper_bounds,
            weights=weights,
            window=model_window,
            get_fit_info=get_fit_info,
            ignore_units=ignore_units,
            **kwargs
        )

        if model_guess.name is not None:
            fit_model.name = model_guess.name

        fitted_models.append(fit_model)

    if single_model_in:
        fitted_models = fitted_models[0]

    return fitted_models


def _fit_lines(spectrum, model, fitter=fitting.LevMarLSQFitter(calc_uncertainties=True),
               exclude_regions=None, exclude_region_upper_bounds=False, weights=None,
               window=None, get_fit_info=False, ignore_units=False, **kwargs):
    """
    Fit the input model (initial conditions) to the spectrum.  Output will be
    the same model with the parameters set based on the fitting.

    spectrum, model -> model
    """
    #
    # If we are to exclude certain regions, then remove them.
    #

    if exclude_regions is not None:
        spectrum = excise_regions(spectrum, exclude_regions,
                                  inclusive_upper=exclude_region_upper_bounds)

    if isinstance(weights, str):
        if weights == 'unc':
            if spectrum.uncertainty is not None:
                # Convert uncertainty to StdDev, then invert
                uncerts = spectrum.uncertainty.represent_as(StdDevUncertainty)
                # Astropy fitters expect weights in 1/sigma
                weights = uncerts.array ** -1
            else:
                weights = None
                warnings.warn("Fitting is set to use uncertainties as weights,"
                              " but the input spectrum's uncertainty is None")
        else:
            raise ValueError("Unrecognized value `%s` in keyword argument.",
                             weights)
    elif weights is not None:
        # Assume that the weights argument is list-like
        weights = np.array(weights)

    mask = spectrum.mask

    dispersion = spectrum.spectral_axis

    flux = spectrum.flux
    flux_unit = spectrum.flux.unit

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
    window_indices = None
    if window is not None and isinstance(window, u.Quantity):
        if isinstance(window.value, (float, int)):
            center = model.mean
            window_indices = np.nonzero((dispersion >= center-window) &
                                        (dispersion < center+window))
        elif len(window.value) == 2:
            window_indices = np.nonzero((dispersion >= window.min()) &
                                        (dispersion <= window.max()))

    # In this case the window is the start and end points of where we
    # should fit
    elif window is not None and isinstance(window, tuple):
        window_indices = np.nonzero((dispersion >= min(window)) &
                                    (dispersion <= max(window)))

    # in this case the window is spectral regions that determine where
    # to fit.
    elif window is not None and isinstance(window, SpectralRegion):
        idx1, idx2 = window.bounds
        if idx1 == idx2:
            raise IndexError("Tried to fit a region containing no pixels.")

        # HACK WARNING! This uses the extract machinery to create a set of
        # indices by making an "index spectrum"
        # note that any unit will do but Jy is at least flux-y
        # TODO: really the spectral region machinery should have the power
        # to create a mask, and we'd just use that...
        idxarr = np.arange(spectrum.flux.size).reshape(spectrum.flux.shape)
        index_spectrum = Spectrum1D(spectral_axis=dispersion,
                                    flux=u.Quantity(idxarr, u.Jy, dtype=int))

        extracted_regions = extract_region(index_spectrum, window)
        if isinstance(extracted_regions, list):
            if len(extracted_regions) == 0:
                raise ValueError('The whole spectrum is windowed out!')
            window_indices = np.concatenate([s.flux.value.astype(int) for s in extracted_regions])
        else:
            if len(extracted_regions.flux) == 0:
                raise ValueError('The whole spectrum is windowed out!')
            window_indices = extracted_regions.flux.value.astype(int)

    if window_indices is not None:
        dispersion = dispersion[window_indices]
        flux = flux[window_indices]
        if mask is not None:
            mask = mask[window_indices]
        if weights is not None:
            weights = weights[window_indices]

    if flux is None or len(flux) == 0:
        raise Exception("Spectrum flux is empty or None.")

    input_spectrum = spectrum

    spectrum = Spectrum1D(
        flux=flux.value * flux_unit,
        spectral_axis=dispersion,
        wcs=input_spectrum.wcs,
        velocity_convention=input_spectrum.velocity_convention,
        rest_value=input_spectrum.rest_value)

    if not model._supports_unit_fitting:
        # Not all astropy models support units.  For those that don't
        # we will strip the units and then re-add them before returning
        # the model.
        model, dispersion, flux = _strip_units_from_model(model, spectrum, convert=not ignore_units)

    #
    # Do the fitting of spectrum to the model.
    #
    if mask is not None:
        nmask = ~mask
        dispersion = dispersion[nmask]
        flux = flux[nmask]
        if weights is not None:
            weights = weights[nmask]

    fit_model = fitter(model, dispersion,
                       flux, weights=weights, **kwargs)

    if hasattr(fitter, 'fit_info') and get_fit_info:
        fit_model.meta['fit_info'] = fitter.fit_info

    if not model._supports_unit_fitting:
        fit_model = QuantityModel(fit_model,
                                  spectrum.spectral_axis.unit,
                                  spectrum.flux.unit)

    return fit_model


def _convert(quantity, dispersion_unit, dispersion, flux_unit):
    """
    Convert the quantity to the spectrum's units, and then we will use
    the *value* of it in the new unitless-model.
    """
    with u.set_enabled_equivalencies(u.spectral()):
        if quantity.unit.is_equivalent(dispersion_unit):
            quantity = quantity.to(dispersion_unit)

    with u.set_enabled_equivalencies(u.spectral_density(dispersion)):
        if quantity.unit.is_equivalent(flux_unit):
            quantity = quantity.to(flux_unit)

    return quantity


def _convert_and_dequantify(poss_quantity, dispersion_unit, dispersion,
                            flux_unit, convert=True):
    """
    This method will convert the ``poss_quantity`` value to the proper
    dispersion or flux units and then strip the units.

    If the ``poss_quantity`` is None, or a number, we just return that.

    Notes
    -----
        This method can be removed along with most of the others here
        when astropy.fitting will fit models that contain units.

    """

    if poss_quantity is None or isinstance(poss_quantity, (float, int)):
        return poss_quantity

    if convert and hasattr(poss_quantity, 'quantity') and poss_quantity.quantity is not None:
        q = poss_quantity.quantity

        quantity = _convert(q, dispersion_unit, dispersion, flux_unit)
        v = quantity.value

    elif convert and isinstance(poss_quantity, u.Quantity):
        quantity = _convert(poss_quantity, dispersion_unit, dispersion, flux_unit)
        v = quantity.value

    else:
        v = poss_quantity.value

    return v


def _strip_units_from_model(model_in, spectrum, convert=True):
    """
    This method strips the units from the model, so the result can
    be passed to the fitting routine. This is necessary as CoumpoundModel
    with units does not work in the fitters.

    Notes
    -----
        When CompoundModel with units works in the fitters this method
        can be removed.

        This assumes there are two types of models, those that are
        based on `~astropy.modeling.models.PolynomialModel` and therefore
        require the ``degree`` parameter when instantiating the class, and
        "everything else" that does not require an "extra" parameter for
        class instantiation.

        If convert is False, then we will *not* do the conversion of units
        to the units of the Spectrum1D object.  Otherwise we will convert.
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
    compound_model = model_in.n_submodels > 1

    if not compound_model:
        # For this we are going to just make it a list so that we
        # can use the looping structure below.
        model_in = [model_in]
    else:
        # If it is a compound model then we are going to create the RPN
        # representation of it which is a list that contains either astropy
        # models or string representations of operators (e.g., '+' or '*').
        model_in = model_in.traverse_postorder(include_operator=True)
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
            new_sub_model = sub_model.__class__(sub_model.degree, name=sub_model.name)
        else:
            new_sub_model = sub_model.__class__(name=sub_model.name)

        # Now for each parameter in the model determine if a dispersion or
        # flux type of unit, then convert to spectrum units and then
        # get the value.

        for pn in new_sub_model.param_names:

            # This could be a Quantity or Parameter
            v = _convert_and_dequantify(getattr(sub_model, pn), dispersion_unit, dispersion, flux_unit, convert=convert)

            #
            # Add this information for the parameter name into the
            # new sub model.
            #
            setattr(new_sub_model, pn, v)

            #
            # Copy over all the constraints (e.g., tied, fixed...)
            #
            for constraint in ('tied', 'fixed'):
                for k, v in getattr(sub_model, constraint).items():
                    getattr(new_sub_model, constraint)[k] = v
            #
            # Convert teh bounds parameter
            #
            new_bounds = []
            for a in sub_model.bounds[pn]:
                v = _convert_and_dequantify(a, dispersion_unit, dispersion, flux_unit, convert=convert)
                new_bounds.append(v)

            new_sub_model.bounds[pn] = tuple(new_bounds)

        # The new model now has unitless information in it but has been
        # converted to spectral unit scale.
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

    Notes
    -----
        When CompoundModel with units works in the fitters this method
        can be removed.

        This assumes there are two types of models, those that are
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

    compound_model = model_in.n_submodels > 1
    if not compound_model:
        model_in_list = [model_in]
        model_orig_list = [model_orig]
    else:
        compound_model_in = model_in

        model_in_list = model_in.traverse_postorder(include_operator=True)
        model_orig_list = model_orig.traverse_postorder(include_operator=True)

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
            new_sub_model = m_in.__class__(m_in.degree, name=m_in.name)
        else:
            new_sub_model = m_in.__class__(name=m_in.name)

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

                    # If it is a compound model, then we need to get the value
                    # from the actual compound model as the tree is not
                    # updated in the fitting
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
                    # If it is a compound model, then we need to get the value
                    # from the actual compound model as the tree is not
                    # updated in the fitting
                    if compound_model:
                        current_value = getattr(compound_model_in, '{}_{}'.format(pn, model_index)).value *\
                                        spectrum.flux.unit
                    else:
                        current_value = m_in_param.value * spectrum.flux.unit

                    v = current_value.to(m_orig_param_quantity.unit,
                                         equivalencies=u.equivalencies.spectral_density(dispersion))
                else:
                    raise ValueError(
                        "The parameter '{}' with unit '{}' is not convertible "
                        "to either the current flux unit '{}' or spectral "
                        "axis unit '{}'.".format(
                            m_orig_param.name, m_orig_param.unit,
                            spectrum.flux.unit, spectrum.spectral_axis.unit))

            else:
                v = getattr(m_in, pn).value

            #
            # Set the parameter value into the new sub-model.
            #

            setattr(new_sub_model, pn, v)

            #
            # Copy over all the constraints (e.g., tied, fixed, bounds...)
            #
            for constraint in ('tied', 'bounds', 'fixed'):
                for k, v in getattr(m_orig, constraint).items():
                    getattr(new_sub_model, constraint)[k] = v

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

    # If the first parameter is not a Quantity, then at this point we will
    # assume none of them are. (It would be inconsistent for fitting to have
    # a model that has some parameters as Quantities and some values).
    if getattr(model_orig, model_orig.param_names[0]).unit is None:
        model_out = QuantityModel(model_out,
                                  spectrum.spectral_axis.unit,
                                  spectrum.flux.unit)

    return model_out


def _combine_postfix(equation):
    """
    Given a Python list in post order (RPN) of an equation, convert/apply the
    operations to evaluate. The list order is the same as what is output from
    ``model._tree.traverse_postorder()``.

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
