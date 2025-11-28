import astropy.units as u
import numpy as np

from astropy.modeling import models, fitting
from specutils.fitting.spline import SplineModel, SplineFitter

from scipy import interpolate


def make_data(with_errs=True):
    """ Arbitrary data """
    np.random.seed(348957)
    x = np.linspace(0, 10, 200)
    y = (x+1) - (x-5)**2. + 10.*np.exp(-0.5 * ((x-7.)/.5)**2.)
    y = (y - np.min(y) + 10.)*10.
    if with_errs:
        ey = np.sqrt(y)
        y = y + np.random.normal(0., ey, y.shape)
    w = 1./y
    return x, y, w


def test_spline_fit():
    x, y, w = make_data()
    make_plot = False

    # Construct three sets of splines and their scipy equivalents
    knots = np.arange(1, 10)
    print(len(x))
    models = [SplineModel(), SplineModel(degree=5), SplineModel(knots=knots),
              SplineModel(smoothing=0)]
    labels = ["Deg 3", "Deg 5", "Knots", "Interpolated"]
    scipyfit = [interpolate.UnivariateSpline(x, y, w),
                interpolate.UnivariateSpline(x, y, w, k=5),
                interpolate.LSQUnivariateSpline(x, y, knots, w=w),
                interpolate.InterpolatedUnivariateSpline(x, y, w)]

    fitter = SplineFitter()
    for model, label, scipymodel in zip(models, labels, scipyfit):
        fitter(model, x, y, w)
        my_y = model(x)
        my_dy = model.derivative()(x)
        my_ady = model.antiderivative()(x)
        my_int = model.integral(x[0],x[-1])
        sci_y = scipymodel(x)
        sci_dy = scipymodel.derivative()(x)
        sci_ady = scipymodel.antiderivative()(x)
        sci_int = scipymodel.integral(x[0],x[-1])
        assert np.allclose(my_y, sci_y, atol=1e-6), label
        assert np.allclose(my_dy, sci_dy, atol=1e-6), label
        assert np.allclose(my_ady, sci_ady, atol=1e-6), label
        assert np.allclose(my_int, sci_int, atol=1e-6), label

        my_ders = model.derivatives(x)
        sci_ders = scipymodel.derivatives(x)
        assert np.allclose(my_ders, sci_ders, atol=1e-6), label
        if model.degree == 3:
            my_roots = model.roots()
            sci_roots = scipymodel.roots()
            assert np.allclose(my_roots, sci_roots, atol=1e-6), label
    
    if make_plot:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.plot(x, y, 'k.')
        ymin, ymax = np.min(y), np.max(y)
        for i, (model, label) in enumerate(zip(models, labels)):
            l, = ax.plot(x, model(x), lw=1, label=label)
            knots = model.knots
            # Hack for now
            if knots is None:
                knots = model._tck[0]

            dy = (ymax-ymin)/10.
            dy /= i+1.
            ax.vlines(knots, ymin, ymin + dy, color=l.get_color(), lw=1)
        ax.legend()
        plt.show()
