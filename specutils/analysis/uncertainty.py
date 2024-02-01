"""
A module for analysis tools dealing with uncertainties or error analysis in
spectra.
"""

import numpy as np
from astropy.nddata import StdDevUncertainty
from ..spectra import SpectralRegion
from ..manipulation import extract_region

__all__ = ['snr', 'snr_derived']


def snr(spectrum, region=None):
    """
    Calculate the mean S/N of the spectrum based on the flux and uncertainty
    in the spectrum. This will be calculated over the regions, if they
    are specified.

    Parameters
    ----------
    spectrum : `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object overwhich the equivalent width will be calculated.

    region: `~specutils.SpectralRegion` or list of `~specutils.SpectralRegion`
        Region within the spectrum to calculate the SNR.

    Returns
    -------
    snr : `~astropy.units.Quantity` or list (based on region input)
        Signal to noise ratio of the spectrum or within the regions

    Notes
    -----
    The spectrum will need to have the uncertainty defined in order for the SNR
    to be calculated. If the goal is instead signal to noise *per pixel*, this
    should be computed directly as ``spectrum.flux / spectrum.uncertainty``. This
    calculation converts the uncertainty to standard deviation internally if it
    is defined as another type.

    """

    if not hasattr(spectrum, 'uncertainty') or spectrum.uncertainty is None:
        raise Exception("Spectrum1D currently requires the uncertainty be defined.")

    # No region, therefore whole spectrum.
    if region is None:
        return _snr_single_region(spectrum)

    # Single region
    elif isinstance(region, SpectralRegion):
        return _snr_single_region(spectrum, region=region)

    # List of regions
    elif isinstance(region, list):
        return [_snr_single_region(spectrum, region=reg)
                for reg in region]


def _snr_single_region(spectrum, region=None):
    """
    Calculate the mean S/N of the spectrum based on the flux and uncertainty
    in the spectrum.

    Parameters
    ----------
    spectrum : `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object overwhich the equivalent width will be calculated.

    region: `~specutils.SpectralRegion`
        Region within the spectrum to calculate the SNR.

    Returns
    -------
    snr : `~astropy.units.Quantity` or list (based on region input)
        Signal to noise ratio of the spectrum or within the regions

    Notes
    -----
    This is a helper function for the above `snr()` method.

    """

    if region is not None:
        calc_spectrum = extract_region(spectrum, region)
    else:
        calc_spectrum = spectrum

    flux = calc_spectrum.flux
    uncertainty = calc_spectrum.uncertainty.represent_as(StdDevUncertainty).quantity

    if hasattr(calc_spectrum, 'mask') and calc_spectrum.mask is not None:
        flux = flux[~calc_spectrum.mask]
        uncertainty = uncertainty[~calc_spectrum.mask]

    # the axis=-1 will enable this to run on single-dispersion, single-flux
    # and single-dispersion, multiple-flux
    return np.mean(flux / uncertainty, axis=-1)


def snr_derived(spectrum, region=None):
    """
    This function computes the signal to noise ratio DER_SNR following the
    definition set forth by the Spectral Container Working Group of ST-ECF,
    MAST and CADC.

    Parameters
    ----------
    spectrum : `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object overwhich the equivalent width will be calculated.

    region: `~specutils.SpectralRegion`
        Region within the spectrum to calculate the SNR.

    Returns
    -------
    snr : `~astropy.units.Quantity` or list (based on region input)
        Signal to noise ratio of the spectrum or within the regions

    Notes
    -----
    The DER_SNR algorithm is an unbiased estimator describing the spectrum
    as a whole as long as the noise is uncorrelated in wavelength bins spaced
    two pixels apart, the noise is Normal distributed, for large wavelength
    regions, the signal over the scale of 5 or more pixels can be approximated
    by a straight line.

    The code and some documentation is derived from
    ``http://www.stecf.org/software/ASTROsoft/DER_SNR/der_snr.py``, and the
    algorithm itself is documented at https://esahubble.org/static/archives/stecfnewsletters/pdf/hst_stecf_0042.pdf
    """

    # No region, therefore whole spectrum.
    if region is None:
        return _snr_derived(spectrum)

    # Single region
    elif isinstance(region, SpectralRegion):
        return _snr_derived(spectrum, region=region)

    # List of regions
    elif isinstance(region, list):
        return [_snr_derived(spectrum, region=reg)
                for reg in region]


def _snr_derived(spectrum, region=None):
    """
    This function computes the signal to noise ratio DER_SNR following the
    definition set forth by the Spectral Container Working Group of ST-ECF,
    MAST and CADC

    Parameters
    ----------
    spectrum : `~specutils.spectra.spectrum1d.Spectrum1D`
        The spectrum object overwhich the equivalent width will be calculated.

    region: `~specutils.SpectralRegion`
        Region within the spectrum to calculate the SNR.

    Returns
    -------
    snr : `~astropy.units.Quantity` or list (based on region input)
        Signal to noise ratio of the spectrum or within the regions

    Notes
    -----
    This is a helper function for the above `snr_derived()` method.

    """

    if region is not None:
        calc_spectrum = extract_region(spectrum, region)
    else:
        calc_spectrum = spectrum

    if hasattr(spectrum, 'mask') and spectrum.mask is not None:
        flux = calc_spectrum.flux[~calc_spectrum.mask]
    else:
        flux = calc_spectrum.flux

    # Values that are exactly zero (padded) are skipped
    n = len(flux)

    # For spectra shorter than this, no value can be returned
    if n > 4:
        signal = np.median(flux)
        noise  = 0.6052697 * np.median(np.abs(2.0 * flux[2:n-2] - flux[0:n-4] - flux[4:n]))
        return signal / noise
    else:
        return 0.0
