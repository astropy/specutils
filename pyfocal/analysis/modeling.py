from astropy.modeling import models, fitting


def apply_model(model, x, y_init, fitter=None):
    if fitter is not None:
        model = fitter(model, x, y_init)

    return model