__author__ = 'Shailesh'

from astropy.modeling import Model, polynomial
from scipy import interpolate

class SplineModel(polynomial.PolynomialModel):
    def __init__(self, degree, knots, coefficients):
        self.knots = knots
        self.coefficients = coefficients
        super(SplineModel, self).__init__(self, degree, n_inputs=1, n_outputs=1, param_dim=1)

    def __call__(self, x):
        return interpolate.splev(x, (self.knots, self.coefficients, self.degree))

