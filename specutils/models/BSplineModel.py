from astropy.modeling import Model, Parameter
from six.moves import xrange

class BSplineModel(Model):
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
        from scipy.interpolate import splrep

        knots, coefficients, _ = splrep(x, y, k=degree)
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
        self.param_names = self._generate_param_names(self.n_pieces)

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
        from scipy.interpolate import splev

        coefficients, knots = [], []
        for i in xrange(self.n_pieces):
            coefficients.append(self.__getattr__("c{:d}".format(i)).value)
            knots.append(self.__getattr__("t{:d}".format(i)).value)
        return splev(x, (knots, coefficients, self.degree))

    def __getattr__(self, attr):
        if self.param_names and attr in self.param_names:
            return Parameter(attr, default=0.0, model=self)
        else:
            return super(BSplineModel, self).__getattr__(attr)

    def __setattr__(self, attr, value):
        # TODO: Support a means of specifying default values for coefficients
        # Check for self._ndim first--if it hasn't been defined then the
        # instance hasn't been initialized yet and self.param_names probably
        # won't work.
        # This has to vaguely duplicate the functionality of
        # Parameter.__set__.
        # TODO: I wonder if there might be a way around that though...
        if attr[0] != '_' and self.param_names and attr in self.param_names:
            param = Parameter(attr, default=0.0, model=self)
            # This is a little hackish, but we can actually reuse the
            # Parameter.__set__ method here
            param.__set__(self, value)
        else:
            super(BSplineModel, self).__setattr__(attr, value)