import numpy as np
import logging
from .resample import Resample
from ..spectra import Spectrum1D


class Splice(object):
    def __init__(self, spacing='coarse'):
        if not spacing in ['coarse', 'fine']:
            raise ValueError("Attribute 'spacing' must be one of 'coarse' or 'fine'.")

        self._spacing = spacing

    def __call__(self, spec_array):
        # Cache resampled Spectrum1D objects
        resamp_specs = []

        # Get the appropriate new bin size for the out dispersion
        fin_disp = None
        min_disp = min([spec.wavelength[0] for spec in spec_array])
        max_disp = max([spec.wavelength[-1] for spec in spec_array])

        for spec in spec_array:
            if self._spacing == 'coarse':
                if fin_disp is None or np.mean(np.diff(spec.wavelength)) > np.mean(np.diff(fin_disp)):
                    fin_disp = spec.wavelength
                    logging.info("Increasing bin width to {}.".format(np.mean(np.diff(spec.wavelength))))
            elif self._spacing == 'fine':
                if fin_disp is None or np.mean(np.diff(spec.wavelength)) < np.mean(np.diff(fin_disp)):
                    fin_disp = spec.wavelength
                    logging.info("Decreasing bin width to {}.".format(np.mean(np.diff(spec.wavelength))))

        fin_disp = np.arange(min_disp.value, max_disp.value, np.mean(np.diff(fin_disp.value))) * fin_disp.unit

        # Define the resample instance
        resample = Resample(fin_disp)

        for spec in spec_array:
            weights = spec.meta.get('weights')

            if weights is None:
                weights = np.ones(spec.flux.shape)

            flux = resample(spec.wavelength, spec.flux)
            uncert = spec.uncertainty.__class__(resample(spec.wavelength,
                                                         spec.uncertainty.array))
            wave = fin_disp
            weights = resample(spec.wavelength, weights)

            resamp_specs.append(Spectrum1D(flux=flux, uncertainty=uncert,
                                           spectral_axis=wave, meta={'weights': weights}))

        return (np.sum([spec * spec.meta['weights'] for spec in resamp_specs], axis=0) /
                np.sum([spec.meta['weights'] for spec in resamp_specs], axis=0))