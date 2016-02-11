from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from astropy.modeling import models, fitting


def apply_model(model, x, y_init, fitter=None):
    if fitter is not None:
        model = fitter(model, x, y_init)

    return model
