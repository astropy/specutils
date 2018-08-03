import numpy as np
from specutils.spectra import Spectrum1D
from astropy.nddata import StdDevUncertainty


def noise_region_uncertainty(spectrum, spectral_region, noise_func=np.std):
    """
    Set the uncertainty in the ``spectrum`` by setting (or overwriting)
    each element of the uncertainty to the the noise of the flux
    calculated from the ``spectral_region``.

    Parameters
    ----------

    spectrum: `~specutils.spectra.Spectrum1D
        The spectrum to which we want to set the uncertainty.

    spectral_region: `~specutils.spectra.SpectralRegion`
        The region to use to calculate the standard deviation.

    Return
    ------
    spectrum_uncertainty: `~specutils.spectra.Spectrum1D
        The spectrum with the uncertainty set.

    """

    # Extract the sub spectrum based on the region
    sub_spectrum = spectral_region.extract(spectrum)

    # Compute the standard deviation of the flux.
    noise = noise_func(sub_spectrum.flux)

    uncertainty = StdDevUncertainty(noise*np.ones(spectrum.flux.shape))

    # Return new specturm with uncertainty set.
    return Spectrum1D(flux=spectrum.flux, spectral_axis=spectrum.spectral_axis,
                      uncertainty=uncertainty,
                      wcs=spectrum.wcs, unit=spectrum.unit,
                      spectral_axis_unit=spectrum.spectral_axis_unit,
                      velocity_convention=spectrum.velocity_convention,
                      rest_value=spectrum.rest_value)
