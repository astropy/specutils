from astropy.nddata import NDDataRef


class NDContainer(NDDataRef):
    """
    Transient n-dimensional data container meant for analysis and manipulation
    functions that require uncertainty propagation without the need to track
    spectral axis and wcs information.
    """
    def __add__(self, other):
        return self.add(other)

    def __sub__(self, other):
        return self.subtract(other)

    def __mul__(self, other):
        return self.multiply(other)

    def __divmod__(self, other):
        return self.divide(other)
