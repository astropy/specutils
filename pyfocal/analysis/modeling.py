from ..interfaces.factories import FitterFactory

def apply_model(model, x, y_init, fitter_name=None):

    if fitter_name:
        fitter = FitterFactory.all_fitters[fitter_name]()
    else:
        fitter = FitterFactory.default_fitter()

    result = fitter(model, x, y_init)

    return result