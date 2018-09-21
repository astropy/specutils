import numpy as np

from astropy import units as u

from .. import Spectrum1D, SpectralRegion
from astropy.nddata import StdDevUncertainty
from .extract_spectral_region import extract_region


__all__ = ['noise_region_uncertainty']


def noise_region_uncertainty(spectrum, spectral_region, noise_func=np.std):
    """
    Generates a new spectrum with an uncertainty from the noise in a particular
    region of the spectrum.

    Parameters
    ----------

    spectrum: `~specutils.spectra.Spectrum1D
        The spectrum to which we want to set the uncertainty.

    spectral_region: `~specutils.spectra.SpectralRegion`
        The region to use to calculate the standard deviation.

    noise_func: callable
        A function which takes the (1D) flux in the ``spectral_region`` and
        yields a *single* value for the noise to use in the result spectrum.

    Return
    ------
    spectrum_uncertainty: `~specutils.spectra.Spectrum1D
        The ``spectrum``, but with a constant uncertainty set by the result of
        the noise region calculation

    """

    # Extract the sub spectrum based on the region
    sub_spectra = extract_region(spectrum, spectral_region)


    # TODO: make this work right for multi-dimensional spectrum1D's?
    if not isinstance(sub_spectra, list):
        sub_spectra = [sub_spectra]
    sub_flux = u.Quantity(np.concatenate([s.flux.value for s in sub_spectra]),
                          spectrum.flux.unit)

    # Compute the standard deviation of the flux.
    noise = noise_func(sub_flux)

    uncertainty = StdDevUncertainty(noise*np.ones(spectrum.flux.shape))

    # Return new specturm with uncertainty set.
    return Spectrum1D(flux=spectrum.flux, spectral_axis=spectrum.spectral_axis,
                      uncertainty=uncertainty,
                      wcs=spectrum.wcs,
                      velocity_convention=spectrum.velocity_convention,
                      rest_value=spectrum.rest_value)
