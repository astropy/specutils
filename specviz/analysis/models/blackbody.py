import numpy as np

import astropy.units as u

from astropy.modeling import Fittable1DModel
from astropy.modeling.parameters import Parameter
from astropy.analytic_functions import blackbody_lambda

__all__ = ['BlackBody']


class BlackBody(Fittable1DModel):

    temp = Parameter(default=5000, min=10.)
    norm = Parameter(default=1.)

    def evaluate(self, x, temp, norm):
        # x is passed as a bare numpy array; must be
        # converted back to Quantity before calling
        # astropy's black body functions.
        _x_u = x * self.wave.unit

        # convert result of the Planck function to
        # flux density in the same units as the data.
        _flux = (blackbody_lambda(_x_u, temp) * u.sr).to(self.flux.unit)

        # normalize and return just the values,
        # to conform to the Model API.
        return (norm * _flux).value


class BlackBodyInitializer(object):

    def initialize(self, instance, wave, flux):
        instance.wave = wave
        instance.flux = flux

        # Roughly normalize by area, with fixed temperature.
        # This may not be good enough for certain types of
        # spectra, so user tweaking may be necessary anyway.
        sum_data = np.sum(flux.value)
        sum_model = np.sum(instance.evaluate(wave.value, instance.temp.value, 1.))

        instance.norm = sum_data / sum_model

