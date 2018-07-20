import astropy.units as u


class SpectralRegion:
    """
    SpectralRegion is currently a container class enables some simplicty 
    in defining and passing a region (interval) for a spectrum.

    In the future, there might be more functionality added in here and there 
    is some discussion that this might/could move to astropy.region.
    """

    def __init__(self, lower, upper):
        """
        Lower and upper values for the interval.

        Parameter
        --------

        lower: number with ``astropy.units`` unit
           The lower bound of the region.

        upper: number with ``astropy.units`` unit
           The upper bound of the region.

        Notes
        -----

        The ``astropy.units`` for each do not have to be the same, but to retain
        sanity it might be easiest to have them the same.  The will need to be
        dispersion units as opposed to flux units.
        """

        if not isinstance(lower, u.Quantity):
            raise TypeError('Lower bound of the region must have an astropy.units unit')

        if not isinstance(upper, u.Quantity):
            raise TypeError('Upper bound of the region must have an astropy.units unit')

        self._lower = lower
        self._upper = upper

    @property
    def lower(self):
        return self._lower

    @lower.setter
    def lower(self, value):
        """
        Lower bound for the interval.

        Parameter
        --------

        value: number with ``astropy.units`` unit
           The lower bound of the region.
        """

        if not isinstance(value, u.Quantity):
            raise TypeError('Lower bound of the region must have an astropy.units unit')

        self._lower = value

    @property
    def upper(self):
        return self._upper

    @upper.setter
    def upper(self, value):
        """
        Upper bound for the interval.

        Parameter
        --------

        value: number with ``astropy.units`` unit
           The upper bound of the region.
        """


        if not isinstance(value, u.Quantity):
            raise TypeError('Upper bound of the region must have an astropy.units unit')

        self._upper = value
