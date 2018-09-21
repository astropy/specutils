from __future__ import division

from astropy.modeling.polynomial import Chebyshev1D
from astropy.modeling.fitting import LevMarLSQFitter

from ..fitting import fit_lines
from ..manipulation.smoothing import median_smooth


__all__ = ['fit_continuum', 'fit_generic_continuum']


def fit_generic_continuum(spectrum, median_window=3, model=Chebyshev1D(3),
                          fitter=LevMarLSQFitter(),
                          exclude_regions=None, weights=None):
    """
    Basic fitting of the continuum of an input spectrum. The input
    spectrum is smoothed using a median filter to remove the spikes.

    Parameters
    ----------
    spectrum : Spectrum1D
        The spectrum object overwhich the equivalent width will be calculated.

    model : list of `~astropy.modeling.Model`
        The list of models that contain the initial guess.

    exclude_regions : list of 2-tuples
        List of regions to exclude in the fitting. Passed through
        to the fitmodels routine.

    weights : list  (NOT IMPLEMENTED YET)
        List of weights to define importance of fitting regions.

    Returns
    -------
    continuum_model
        Fitted continuum as a model of whatever class ``model`` provides.

    Notes
    -----
       * Could add functionality to set the bounds in
         ``model`` if they are not set.
       * The models in the list of ``model`` are added together and passed as a
         compound model to the `~astropy.modeling.fitting.Fitter` class instance.

    """

    #
    # Simple median smooth to remove spikes and peaks
    #

    spectrum_smoothed = median_smooth(spectrum, median_window)

    #
    # Return the fitted continuum
    #

    return fit_continuum(spectrum_smoothed, model, fitter, exclude_regions, weights)


def fit_continuum(spectrum, model=Chebyshev1D(3), fitter=LevMarLSQFitter(),
                  exclude_regions=None, window=None, weights=None):
    """
    Entry point for fitting using the `~astropy.modeling.fitting`
    machinery.

    Parameters
    ----------
    spectrum : Spectrum1D
        The spectrum object overwhich the equivalent width will be calculated.

    model: list of `~astropy.modeling.Model`
        The list of models that contain the initial guess.

    fitmodels_type: str
        String representation of fit method to use as defined by the dict fitmodels_types.

    window : tuple of wavelengths  (NOT IMPLEMENTED YET)
        Start and end wavelengths used for fitting.

    weights : list  (NOT IMPLEMENTED YET)
        List of weights to define importance of fitting regions.

    Returns
    -------
    models : list of `~astropy.modeling.Model`
        The list of models that contain the fitted model parmeters.

    """

    if window is not None or weights is not None:
        raise NotImplementedError('window and weights are not yet implemented')

    #
    # Fit the flux to the model.
    #

    continuum_spectrum = fit_lines(spectrum, model, fitter, exclude_regions, weights)

    return continuum_spectrum
