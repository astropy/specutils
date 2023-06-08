"""
A module for analysis tools focused on determining the location of
spectral features.
"""
import warnings

from astropy.nddata import StdDevUncertainty
from astropy.utils.exceptions import AstropyDeprecationWarning
import astropy.uncertainty as unc
import numpy as np

from ..spectra import SpectralRegion
from ..manipulation import extract_region


__all__ = ['centroid']


def centroid(spectrum, regions=None, region=None, analytic=True):
    """
    Calculate the centroid of a region, or regions, of the spectrum.

    Parameters
    ----------
    spectrum : `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object over which the centroid will be calculated.  If the uncertainty
        is populated, the returned quantity will include an uncertainty attribute with
        the propagated uncertainty (as Standard Deviation-style uncertainties).  This uncertainty
        assumes the input uncertainties are uncorrelated.

    regions : `~specutils.SpectralRegion` or list of `~specutils.SpectralRegion`
        Region within the spectrum to calculate the centroid.

    analytic : bool, optional
        Set this flag to ``False`` to use the `~astropy.uncertainty` distribution-based
        calculation for the centroid and its uncertainty instead of the default analytic solution.

    Returns
    -------
    centroid : float or list (based on region input)
        Centroid of the spectrum or within the regions

    Notes
    -----
    The spectrum will need to be continuum subtracted before calling
    this method. See the
    `analysis documentation <https://specutils.readthedocs.io/en/latest/analysis.html>`_ for more information.

    """

    if region is not None:
        regions = region
        warnings.warn("The 'region' keyword has been deprecated in favor "
                      "of 'regions' since specutils 1.8 and will be removed "
                      "in a future release.", AstropyDeprecationWarning)

    # No region, therefore whole spectrum.
    if regions is None:
        return _centroid_single_region(spectrum, analytic=analytic)

    # Single region
    elif isinstance(regions, SpectralRegion):
        return _centroid_single_region(spectrum, region=regions,
                                       analytic=analytic)

    # List of regions
    elif isinstance(regions, list):
        return [_centroid_single_region(spectrum, region=reg,
                                        analytic=analytic)
                for reg in regions]


def _centroid_single_region(spectrum, region=None, analytic=True):
    """
    Calculate the centroid of the spectrum based on the flux in the spectrum.
    The returned quantity object will have a ``.uncertainty`` attribute which
    will be populated if ``spectrum`` has uncertainties assigned, or ``None`` if not.

    Parameters
    ----------
    spectrum : `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object overwhich the centroid will be calculated.

    region : `~specutils.SpectralRegion`
        Region within the spectrum to calculate the centroid.

    analytic : bool, optional
        Set this flag to ``False`` to use the `~astropy.uncertainty` distribution-based
        calculation for the centroid and its uncertainty instead of the default analytic solution.

    Returns
    -------
    centroid : float or list (based on region input)
        Centroid of the spectrum or within the regions

    Notes
    -----
    This is a helper function for the above `centroid()` method.

    """
    if not analytic and spectrum.uncertainty is None:
        raise ValueError("Distribution-based calculation can only be used if"
                         " spectrum.uncertainty is not None")

    if region is not None:
        calc_spectrum = extract_region(spectrum, region)
    else:
        calc_spectrum = spectrum

    if spectrum.uncertainty is not None:
        flux_uncert = calc_spectrum.uncertainty.represent_as(StdDevUncertainty).quantity
    else:
        # dummy value for uncertainties to avoid extra if-statements when applying mask
        flux_uncert = np.zeros_like(calc_spectrum.flux)

    if hasattr(calc_spectrum, 'mask') and calc_spectrum.mask is not None:
        flux = calc_spectrum.flux[~calc_spectrum.mask]
        dispersion = calc_spectrum.spectral_axis[~calc_spectrum.mask].quantity
        flux_uncert = flux_uncert[~calc_spectrum.mask]
    else:
        flux = calc_spectrum.flux
        dispersion = calc_spectrum.spectral_axis.quantity

    if analytic:
        centroid = np.sum(flux*dispersion, axis=-1) / np.sum(flux, axis=-1)
        if spectrum.uncertainty is None:
            centroid.uncertainty = None
        else:
            N = np.sum(flux, axis=-1)
            # Looks overcomplicated, but gives us the right shape to match flux_uncert
            diff = np.subtract.outer(dispersion, centroid.transpose()).transpose()
            s2 = np.sum(flux_uncert*flux_uncert * diff*diff, axis=-1)/(N*N)

            # centroid uncertainty, fractionally added in quadrature of numerator and denom
            centroid.uncertainty = np.sqrt(s2).to(spectrum.spectral_axis.unit)
            centroid.uncertainty_type = 'stddev'
    else:
        # Convert flux to an astropy.uncertainties normal distribution
        flux = unc.normal(flux, std=flux_uncert, n_samples=1000)

        if len(flux.shape) > 1:
            dispersion = np.tile(dispersion, list(flux.shape[:-1]) + [1,])

        # centroid is the flux-weighted mean of the dispersions, the uncertainties
        # need to be scaled to the numerator/denominator, so we'll compute those in advance.

        # the axis=-2 will enable this to run on single-dispersion, single-flux
        # and single-dispersion, multiple-flux. Not -1 because of the distribution
        # n_samples axis
        numerator = (flux*dispersion).sum(axis=-2)
        denom = flux.sum(axis=-2)
        centroid_dist = numerator/denom
        centroid = centroid_dist.pdf_mean()
        if spectrum.uncertainty is None:
            centroid.uncertainty = None
        else:
            centroid.uncertainty = centroid_dist.pdf_std()
            centroid.uncertainty_type = 'stddev'

    return centroid
