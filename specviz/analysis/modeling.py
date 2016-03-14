from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import logging

from ..interfaces import factories as factories


def apply_model(model, x, y_init, fitter_name=None):

    if fitter_name:
        fitter = factories.FitterFactory.all_fitters[fitter_name]()
    else:
        fitter = factories.FitterFactory.default_fitter()

    result = fitter(model, x, y_init)

    if 'message' in fitter.fit_info:
        # The fitter 'message' should probably be logged at INFO level.
        # Problem is, info messages do not display in the error console,
        # and we, ideally, want the user to see the message immediately
        # after the fit is executed.
        logging.warning(fitter.fit_info['message'])

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
