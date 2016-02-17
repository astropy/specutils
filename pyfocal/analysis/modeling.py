from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ..interfaces.factories import FitterFactory


def apply_model(model, x, y_init, fitter_name=None):

    if fitter_name:
        fitter = FitterFactory.all_fitters[fitter_name]()
    else:
        fitter = FitterFactory.default_fitter()

    result = fitter(model, x, y_init)

    return result


# def gaussian(x, y):
#     amp, mean, stddev = _gaussian_parameter_estimates(x, y)
#     g_init = models.Gaussian1D(amplitude=amp, mean=mean, stddev=stddev)
#     fit_g = fitting.LevMarLSQFitter()
#     g = fit_g(g_init, x, y)
#
#     return (g.amplitude, g.mean, g.stddev), x, g(x)
#
#
# def _gaussian_parameter_estimates(x, y, dy=0):
#     amplitude = np.percentile(y, 95)
#     y = np.max(y / y.sum(), 0)
#     mean = (x * y).sum()
#     stddev = np.sqrt((y * (x - mean) ** 2).sum())
#     return amplitude, mean, stddev
