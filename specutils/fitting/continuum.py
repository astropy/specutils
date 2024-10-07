import astropy.units as u
from astropy.modeling.polynomial import Chebyshev1D
from astropy.modeling.fitting import TRFLSQFitter

from ..fitting import fit_lines
from ..manipulation.smoothing import median_smooth
from ..spectra import SpectralRegion


__all__ = ['fit_continuum', 'fit_generic_continuum']


def fit_generic_continuum(spectrum, median_window=3, model=Chebyshev1D(3),
                          fitter=TRFLSQFitter(),
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
    median_window : float
        The width of the median smoothing kernel used to filter the data before
        fitting the continuum. See the ``kernel_size`` parameter of
        `~scipy.signal.medfilt` for more information.
    fitter : `~astropy.fitting._FitterMeta`
        The astropy fitter to use for fitting the model.
        Default: `~astropy.modeling.fitting.TRFLSQFitter`
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
    # Simple median smooth to remove spikes and peaks
    spectrum_smoothed = median_smooth(spectrum, median_window)

    return fit_continuum(spectrum_smoothed, model=model, fitter=fitter,
                         exclude_regions=exclude_regions, weights=weights)


def fit_continuum(spectrum, model=Chebyshev1D(3), fitter=TRFLSQFitter(),
                  exclude_regions=None, exclude_region_upper_bounds=False,
                  window=None, weights=None):
    """
    Entry point for fitting using the `~astropy.modeling.fitting`
    machinery.

    Parameters
    ----------
    spectrum : Spectrum1D
        The spectrum object overwhich the equivalent width will be calculated.
    model: list of `~astropy.modeling.Model`
        The list of models that contain the initial guess.
    fitter : `~astropy.fitting._FitterMeta`
        The astropy fitter to use for fitting the model.
        Default: `~astropy.modeling.fitting.TRFLSQFitter`
    exclude_regions : list of 2-tuples
        List of regions to exclude in the fitting. Passed through
        to the fitmodels routine.
    exclude_region_upper_bounds : bool, optional
        Set to True to override the default pythonic slicing on the excluded regions (inclusive
        lower bound, exclusive upper bound) and remove spectral values less than or equal to
        the upper bound of the region rather than strictly less than the upper bound. Passed
        through to the fitmodels routine.
    window : tuple of wavelengths
        Start and end wavelengths used for fitting.
    weights : list  (NOT IMPLEMENTED YET)
        List of weights to define importance of fitting regions.

    Returns
    -------
    models : list of `~astropy.modeling.Model`
        The list of models that contain the fitted model parameters.
    """
    if weights is not None:
        raise NotImplementedError("Weights support is not yet implemented.")

    w = window
    if type(w) in [list, tuple] and all([_is_valid_sequence(x) for x in w]):
        w = SpectralRegion(w)

    # Fit the flux to the model
    continuum_spectrum = fit_lines(spectrum, model, fitter, exclude_regions,
                                   exclude_region_upper_bounds, weights, w)

    return continuum_spectrum


# Checks for sequences of of 2-tuples with Quantities
def _is_valid_sequence(value):
    if type(value) in [list, tuple]:
        return len(value) == 2 and \
               isinstance(value[0], u.Quantity) and \
               isinstance(value[1], u.Quantity)

    return False
