from astropy.modeling import polynomial, Parameter
from scipy import interpolate
import numpy as np

class BSplineModel(polynomial.PolynomialBase):
    """
    Implements a BSpline representation of 1-D curve.

    Parameters
    -----------
    degree : int
        the degree of the spline (e.g. 1 for linear, 3 for cubic)
    knots : array_like
        array of knots of the spline
    coefficients: array_like
        coefficients corresponding to the knots, should be
        the same length as knots

    Raises
    --------
    ValueError
        If the length of knots and coefficients arrays do not match
    """

    @classmethod
    def from_data(cls, x, y, degree):
        """
        Initializes the B-spline representation of 1-D curve.
        Given the set of data points (x[i], y[i]) determines a smooth spline
        approximation
        Parameters
        -----------
        x, y : array_like
            The data points defining a curve y = f(x)
        degree: int
            the degree of the spline (e.g. 1 for linear, 3 for cubic)
        """
        knots, coefficients, _ = interpolate.splrep(x, y, k=degree)
        return cls(degree, knots, coefficients)

    def __init__(self, degree, knots, coefficients):
        if len(knots) != len(coefficients):
            raise ValueError("Number of knots ({0}) have to be equal to the"
                             "number of coefficients ({1})"
                             .format(len(knots), len(coefficients)))
        # TODO: Scipy does not raise any error in any case, figure out why
        # if len(knots) - degree - 1 <= 1:
        #     raise ValueError("Not enough knots ({0})".format(len(knots)))

        self.degree = degree
        self.n_pieces = len(knots)
        self._param_names = self._generate_param_names(self.n_pieces)

        params = {}
        for i in xrange(len(knots)):
            params["c{:d}".format(i)] = coefficients[i]
            params["t{:d}".format(i)] = knots[i]
        super(BSplineModel, self).__init__(param_dim=1, **params)

    def _generate_param_names(self, n_pieces):
        names = []
        for i in xrange(n_pieces):
            names.append("c{:d}".format(i))
            names.append("t{:d}".format(i))
        return names

    def __call__(self, x):
        coefficients, knots = [], []
        for i in xrange(self.n_pieces):
            coefficients.append(self.__getattr__("c{:d}".format(i)).value)
            knots.append(self.__getattr__("t{:d}".format(i)).value)
        return interpolate.splev(x, (knots, coefficients, self.degree))

# TODO: convert to a test
if __name__ == "__main__":
    knots = [0., 0., 0., 0., 2.22222222, 3.33333333, 4.44444444, 5.55555556,
             6.66666667, 7.77777778, 10., 10., 10., 10.]
    coefficients = [ -4.94881722e-18,   8.96543619e-01,   1.39407154e+00,
                    -2.36640266e-01,  -1.18324030e+00,  -8.16301228e-01,
                    4.57836125e-01,   1.48720677e+00,   1.64338775e-01,
                    -5.44021111e-01,   0.00000000e+00,   0.00000000e+00,
                    0.00000000e+00,   0.00000000e+00]
    degree = 3
    bs = BSplineModel(degree, knots, coefficients)
    x = np.linspace(0, 10, 100)
    print bs(x)