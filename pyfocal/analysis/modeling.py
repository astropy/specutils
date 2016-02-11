from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from astropy.modeling import models, fitting
from ..interfaces.factories import FitterFactory

def apply_model(model, x, y_init, fitter_name=None):

    if fitter_name:
        fitter = FitterFactory.all_fitters[fitter_name]()
    else:
        fitter = FitterFactory.default_fitter()

    return model