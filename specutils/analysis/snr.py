import numpy as np
from ..spectra import SpectralRegion

__all__ = ['snr']


def snr(spectrum, region=None, noise_region=None):
    """
    Calculate the mean S/N of the spectrum based on the flux and uncertainty in the spectrum. This will
    be calculated over the regions, if they are specified.

    Parameters
    ----------
    spectrum : ``specutils.spectra.spectrum1d.Spectrum1D``
        The spectrum object overwhich the equivalent width will be calculated.

    region: ``specutils.utils.SpectralRegion`` or list of ``specutils.utils.SpectralRegion``
        Region within the spectrum to calculate the SNR.

    noise_region: ``specutils.utils.SpectralRegion``
        Region within the spectrum from which to calculate the noise (based on standard deviation
        of the flux values)

    Returns
    -------
    snr : float or list (based on region input)
        Signal to noise ratio of the spectrum or within the regions

    Notes
    -----
    The spectrum will need to have the uncertainty defined in order
    for the SNR to be calculated.

    Two methods implemented:
      1) SNR based on uncertainty
          This will compute the mean of the flux / uncertainty.

      2) SNR based on noise_region
          This will compute the standard devation of the flux in the noise_region
          and the mean of the flux in the region and then return the mean flux divided
          by the standard deviation of the noise_region flux.

    """

    # No region, therefore whole spectrum.
    if region is None:
        return _snr_single_region(spectrum, noise_region=noise_region)

    # Single region
    elif isinstance(region, SpectralRegion):
        return _snr_single_region(spectrum, region=region, noise_region=noise_region)

    # List of regions
    elif isinstance(region, list):
        return [_snr_single_region(spectrum, region=reg, noise_region=noise_region)
                for reg in region]


def _snr_single_region(spectrum, region=None, noise_region=None):
    """
    Calculate the mean S/N of the spectrum based on the flux and uncertainty in the spectrum.

    Parameters
    ----------
    spectrum : ``specutils.spectra.spectrum1d.Spectrum1D``
        The spectrum object overwhich the equivalent width will be calculated.

    region: ``specutils.utils.SpectralRegion``
        Region within the spectrum to calculate the SNR.

    noise_region: ``specutils.utils.SpectralRegion``
        Region within the spectrum from which to calculate the noise (based on standard deviation
        of the flux values)

    Returns
    -------
    snr : float or list (based on region input)
        Signal to noise ratio of the spectrum or within the regions

    Notes
    -----
    This is a helper function for the above `snr()` method.

    """

    if region is not None:
        calc_spectrum = region.extract(spectrum)
    else:
        calc_spectrum = spectrum

    flux = calc_spectrum.flux

    # If a noise_region is defined then calculate the standard deviation
    # in the noise_region and use that as the uncertainty.
    if noise_region is not None:
        spectrum_noise_region = noise_region.extract(spectrum)
        uncertainty = np.std(spectrum_noise_region.flux)

    # If noise_region is not defined then use the Uncertainty in the spectrum.
    else:
        uncertainty = calc_spectrum.uncertainty.array * spectrum.uncertainty.unit

    return np.mean(flux / uncertainty)
