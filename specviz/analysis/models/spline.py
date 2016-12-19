from scipy.interpolate import UnivariateSpline

from astropy.modeling import Fittable1DModel
from astropy.modeling.parameters import Parameter

__all__ = ['Spline1D']


class Spline1D(Fittable1DModel):
    """
    Perform a spline fit

    Notes
    -----
    See
    `astropy.modeling.Fittable1DModel <http://docs.astropy.org/en/stable/api/astropy.modeling.Fittable1DModel.html#astropy.modeling.Fittable1DModel>`_
    for further details.

    The spline function is based on
    `scipy.interpolate.UnivariateSpline <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.UnivariateSpline.html>`_
    """
    degree = Parameter(default=3, fixed=True)
    smooth = Parameter(default=1, fixed=True)  # default=None crashes the app

    def evaluate(self, x, degree, smooth):
        """
        Evaluate the spline

        Parameters
        ----------
        x: numpy.ndarray
            The wavelengths to evaluate over.

        degree: int
            The degree of spline to evaluate.

        smooth: float or None
            The smoothing factor used to choose the number of knots.

        Returns
        -------
        The evaluated spline
        """
        _f = UnivariateSpline(self.wave, self.flux,
                              k=degree,
                              s=smooth)
        return _f(x)


class Spline1DInitializer(object):
    """
    `Spline1D` model initializer
    """

    def initialize(self, instance, wave, flux):
        """
        Initialize the `Spline1D` model.

        Parameters
        ----------
        instance: `Spline1D` instance
            The `Spline1D` model.

        wave: numpy.ndarray
            The wavelength range.

        flux: numpy.ndarray
            The source flux.
        """
        instance.wave = wave
        instance.flux = flux

        # these override the defaults to something sensible.
        instance.smooth.value = len(wave)
