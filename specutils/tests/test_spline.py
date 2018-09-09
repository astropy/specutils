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
    make_plot=False
    
    # Construct three sets of splines and their scipy equivalents
    knots = np.arange(1,10)
    models = [SplineModel(), SplineModel(degree=5), SplineModel(knots=knots), SplineModel(smoothing=0)]
    labels = ["Deg 3", "Deg 5", "Knots", "Interpolated"]
    scipyfit = [interpolate.UnivariateSpline(x,y,w),
                interpolate.UnivariateSpline(x,y,w,k=5),
                interpolate.LSQUnivariateSpline(x,y,knots,w=w),
                interpolate.InterpolatedUnivariateSpline(x,y,w)]
    
    fitter = SplineFitter()
    for model, label, scipymodel in zip(models, labels, scipyfit):
        fitter(model, x, y, w)
        my_y = model(x)
        sci_y = scipymodel(x)
        assert np.allclose(my_y, sci_y, atol=1e-6)
    
    if make_plot:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.plot(x,y,'k.')
        ymin, ymax = np.min(y), np.max(y)
        for i,(model, label) in enumerate(zip(models, labels)):
            l, = ax.plot(x, model(x), lw=1, label=label)
            knots = model.get_knots()
            # Hack for now
            if knots is None: knots = model._tck[0]
            print(knots)
            dy = (ymax-ymin)/10.
            dy /= i+1.
            ax.vlines(knots, ymin, ymin + dy, color=l.get_color(), lw=1)
        ax.legend()
        plt.show()

if __name__=="__main__":
    test_spline_fit()

