from scipy.interpolate import UnivariateSpline

from astropy.modeling import Fittable1DModel
from astropy.modeling.parameters import Parameter

__all__ = ['Spline1D']


class Spline1D(Fittable1DModel):

    degree = Parameter(default=3, fixed=True)
    smooth = Parameter(default=1, fixed=True)  # default=None crashes the app

    def evaluate(self, x, degree, smooth):
        _f = UnivariateSpline(self.wave, self.flux,
                              k=degree,
                              s=smooth)
        return _f(x)


class Spline1DInitializer(object):
    
    def initialize(self, instance, wave, flux):
        instance.wave = wave
        instance.flux = flux

        # these override the defaults to something sensible.
        instance.smooth.value = len(wave)
