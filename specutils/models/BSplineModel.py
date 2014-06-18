from astropy.modeling import Model, Parameter
from six.moves import xrange

class BSplineModel(Model):
    """
    Implements a BSpline representation of 1-D curve.

    Parameters
    -----------
    degree : int
        the degree of the spline (e.g. 1 for linear, 3 for cubic)
    x : array_like
        If fitted is False, then x defines the data points defining a curve y
        = f(x). Otherwise x is an array of knots of the spline
    y : array_like
        If fitted is False, then y defines the data points defining a curve y
        = f(x). Otherwise y is an array of coefficients corresponding to the
        knots, should be the same length as knots
    fitted : boolean
        If fitted is True, then x and y are treated as knots and coefficients
        respectively. Otherwise, a smooth spline approximation is determined
        from the given set of data points (x[i], y[i])
    Raises
    --------
    ValueError
        If the length of x (knots) and y (coefficients) arrays do not match
    """

    def __init__(self, degree, x, y, fitted=False):

        from scipy.interpolate import splrep
        if len(x) != len(y):
            raise ValueError("Number of elements in x (knots) ({0}) have to be "
                             "equal to the number of elements in y "
                             "(coefficients) ({1})".format(len(x), len(y)))
        # TODO: Scipy does not raise any error in any case, figure out why
        # if len(knots) - degree - 1 <= 1:
        #     raise ValueError("Not enough knots ({0})".format(len(knots)))

        if fitted:
            knots, coefficients = x, y
        else:
            knots, coefficients, _ = splrep(x, y, k=degree)

        self.degree = degree
        self.length = len(knots)
        self.param_names = self._generate_param_names()

        params = {}
        for i in xrange(len(knots)):
            params["c{:d}".format(i)] = coefficients[i]
            params["t{:d}".format(i)] = knots[i]
        super(BSplineModel, self).__init__(param_dim=1, **params)

    def _generate_param_names(self):
        names = []
        for i in xrange(self.length):
            names.append("c{:d}".format(i))
            names.append("t{:d}".format(i))
        return names

    def __call__(self, x):
        from scipy.interpolate import splev

        coefficients, knots = [], []
        for i in xrange(self.length):
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