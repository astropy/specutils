import itertools
import operator
import logging

import numpy as np
from astropy.modeling import fitting, Model, models
from astropy.table import QTable
from scipy.signal import convolve


import astropy.units as u

from ..spectra.spectral_region import SpectralRegion
from ..spectra.spectrum1d import Spectrum1D
from ..analysis import fwhm, gaussian_sigma_width, centroid, warn_continuum_below_threshold
from ..manipulation import extract_region, noise_region_uncertainty
from ..manipulation.utils import excise_regions

__all__ = ['find_lines_threshold', 'find_lines_derivative', 'fit_lines',
           'estimate_line_parameters']

log = logging.getLogger(__name__)


# Define the initial estimators. This are the default methods to use to
# estimate astropy model parameters. This is based on only a small subset of
# the astropy models but it was determined that this is a decent start as most
# fitting will probably use one of these.
#
# Each method list must take a Spectrum1D object and should return a Quantity.
_parameter_estimators = {
    'Gaussian1D': {
        'amplitude': lambda s: max(s.flux),
        'mean': lambda s: centroid(s, region=None),
        'stddev': lambda s: gaussian_sigma_width(s)
    },
    'Lorentz1D': {
        'amplitude': lambda s: max(s.flux),
        'x_0': lambda s: centroid(s, region=None),
        'fwhm': lambda s: fwhm(s)
    },
    'Voigt1D': {
        'x_0': lambda s: centroid(s, region=None),
        'amplitude_L': lambda s: max(s.flux),
        'fwhm_L': lambda s: fwhm(s) / np.sqrt(2),
        'fwhm_G': lambda s: fwhm(s) / np.sqrt(2)
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


def estimate_line_parameters(spectrum, model):
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
            setattr(model, name, estimator(spectrum))
        except AttributeError:
            raise Exception('No method to estimate parameter {}'.format(name))

    return model


def _consecutive(data, stepsize=1):
    return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)


@warn_continuum_below_threshold(threshold=0.01)
def find_lines_threshold(spectrum, noise_factor=1):
    """
    Find the emission and absorption lines in a spectrum. The method
    here is based on deviations larger than the spectrum's uncertainty by the
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
    uncert_val = 0 if uncertainty is None else uncertainty.array

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
    spectrum : Spectrum1D
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


def fit_lines(spectrum, model, fitter=fitting.LevMarLSQFitter(),
              exclude_regions=None, weights=None, window=None,
              **kwargs):
    """
    Fit the input models to the spectrum. The parameter values of the
    input models will be used as the initial conditions for the fit.

    Parameters
    ----------
    spectrum : Spectrum1D
        The spectrum object over which the equivalent width will be calculated.
    model: `~astropy.modeling.Model` or list of `~astropy.modeling.Model`
        The model or list of models that contain the initial guess.
    fitter : `~astropy.modeling.fitting.Fitter`, optional
        Fitter instance to be used when fitting model to spectrum.
    exclude_regions : list of `~specutils.SpectralRegion`
        List of regions to exclude in the fitting.
    weights : array-like or 'unc', optional
        If 'unc', the unceratinties from the spectrum object are used to
        to calculate the weights. If array-like, represents the weights to
        use in the fitting.  Note that if a mask is present on the spectrum, it
        will be applied to the ``weights`` as it would be to the spectrum
        itself.
    window : `~specutils.SpectralRegion` or list of `~specutils.SpectralRegion`
        Regions of the spectrum to use in the fitting. If None, then the
        whole spectrum will be used in the fitting.
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
        spectrum = excise_regions(spectrum, exclude_regions)

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

        fit_model = _fit_lines(spectrum, model_guess, fitter,
                               exclude_regions, weights, model_window,
                               ignore_units, **kwargs)
        if model_guess.name is not None:
            fit_model.name = model_guess.name

        fitted_models.append(fit_model)

    if single_model_in:
        fitted_models = fitted_models[0]

    return fitted_models


def _fit_lines(spectrum, model, fitter=fitting.LevMarLSQFitter(),
               exclude_regions=None, weights=None, window=None,
               ignore_units=False, **kwargs):
    """
    Fit the input model (initial conditions) to the spectrum.  Output will be
    the same model with the parameters set based on the fitting.

    spectrum, model -> model
    """
    #
    # If we are to exclude certain regions, then remove them.
    #

    if exclude_regions is not None:
        spectrum = excise_regions(spectrum, exclude_regions)

    if isinstance(weights, str):
        if weights == 'unc':
            uncerts = spectrum.uncertainty

            # Astropy fitters expect weights in 1/sigma
            if uncerts is not None:
                weights = uncerts.array ** -1
            else:
                log.warning("Uncertainty values are not defined, but are "
                                "trying to be used in model fitting.")
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
    if window is not None and isinstance(window, (float, int)):
        center = model.mean
        window_indices = np.nonzero((dispersion >= center-window) &
                                    (dispersion < center+window))

    # In this case the window is the start and end points of where we
    # should fit
    elif window is not None and isinstance(window, tuple):
        window_indices = np.nonzero((dispersion >= window[0]) &
                                    (dispersion <= window[1]))

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

    #
    # Do the fitting of spectrum to the model.
    #
    if mask is not None:
        nmask = ~mask
        dispersion = dispersion[nmask]
        flux = flux[nmask]
        if weights is not None:
            weights = weights[nmask]

    return fitter(model, dispersion,
                  flux, weights=weights, **kwargs)

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
