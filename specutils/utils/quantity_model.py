from astropy import units as u

__all__ = ['QuantityModel']


class QuantityModel:
    """
    The QuantityModel was created to wrap `~astropy.modeling.models` that do not have
    the ability to use `~astropy.units` in the parameters.

    Parameters
    ----------
    unitless_model: `~astropy.modeling.Model`
       A model that does not have units

    input_units: `~astropy.units`
      Units for the dispersion axis

    return_units: `~astropy.units`
      Units for the flux axis

    Notes
    -----
    When Astropy's modeling is updated so *all* models have the ability
    to have `~astropy.units.Quantity` on all parameters, then this will
    not be needed.
    """

    def __init__(self, unitless_model, input_units, return_units):
        self.unitless_model = unitless_model  # should  check that it's unitless somehow!

        # we use the dict because now this "shadows" the unitless model's input_units/ return_units
        self.__dict__['input_units'] = input_units
        self.__dict__['return_units'] = return_units

    def __hasattr_(self, nm):
        if nm in self.__dict__:
            return True
        if hasattr(self, self.unitless_model):
            return True
        return False

    def __getattr__(self, nm):
        if hasattr(self.unitless_model, nm):
            return getattr(self.unitless_model, nm)
        else:
            raise AttributeError("'{}' object has no attribute '{}'".format(self.__class__.__name__, nm))

    def __setattr__(self, nm, val):
        if nm != 'unitless_model' and hasattr(self.unitless_model, nm):
            setattr(self.unitless_model, nm, val)
        else:
            super().__setattr__(nm, val)

    def __delattr__(self, nm):
        if hasattr(self.unitless_model, nm):
            delattr(self.unitless_model, nm, val)
        else:
            super().__delattr__(nm, val)

    def __dir__(self):
        thisdir = super().__dir__()
        modeldir = dir(self.unitless_model)
        return sorted(thisdir + modeldir)

    def __repr__(self):
        return ('<QuantityModel {}, input_units={}, '
                'return_units={}>'.format(repr(self.unitless_model)[1:-1],
                                          self.input_units, self.return_units))

    def __call__(self, x, *args, **kwargs):
        unitlessx = x.to(self.input_units).value
        result = self.unitless_model(unitlessx, *args, **kwargs)
        return u.Quantity(result, self.return_units, copy=False)
