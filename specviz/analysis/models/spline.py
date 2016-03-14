from scipy.interpolate import UnivariateSpline

from astropy.modeling import Fittable1DModel
from astropy.modeling.parameters import Parameter

__all__ = ['Spline1D']


class Spline1D(Fittable1DModel):

    degree = Parameter(default=3, min=0, max=5, fixed=True)

    def evaluate(self, x, degree):
        return UnivariateSpline(self.wave, self.flux, k=degree, s=100000)(x)


class Spline1DInitializer(object):

    def initialize(self, instance, wave, flux):
        instance.wave = wave
        instance.flux = flux
