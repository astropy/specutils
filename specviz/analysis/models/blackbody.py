import numpy as np

from astropy.modeling import Fittable1DModel
from astropy.modeling.parameters import Parameter
from astropy.analytic_functions import blackbody_lambda, blackbody_nu

__all__ = ['BlackBody']


class BlackBody(Fittable1DModel):

    temp = Parameter(default=5000, min=10.)
    norm = Parameter(default=1.)

    def evaluate(self, x, temp, norm):
        # x is passed as a bare numpy array; must be
        # converted back to Quantity before calling
        # astropy's black body functions.
        _x_u = x * self.wave.unit

        # call the Planck function most appropriate
        # for the flux units being used. Is there a
        # better way to tell apart flam from fnu?
        if str(self.flux.unit).lower().index("angstrom") > 0:
            _flux = blackbody_lambda(_x_u, temp)
        else:
            _flux = blackbody_nu(_x_u, temp)

        # In practice it's not necessary to convert to
        # flux density (by removing the /sr part from the
        # unit) since we normalize to the data at hand
        # anyway. And don't forget to return just the
        # values so as to conform to the Model API.
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

