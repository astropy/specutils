import astropy.units as u


class SpectralRegion:

    def __init__(self, lower, upper):

        self._lower = lower
        self._upper = upper

    @property
    def lower(self):
        return self._lower

    @lower.setter
    def lower(self, value):

        if not isinstance(value, u.Quantity):
            raise TypeError('Lower bound of the region must have an astropy.units unit')

        self._lower = value

    @property
    def upper(self):
        return self._upper

    @upper.setter
    def upper(self, value):

        if not isinstance(value, u.Quantity):
            raise TypeError('Upper bound of the region must have an astropy.units unit')

        self._upper = value
