from specuctils import Spectrum1D


class RadioSpectrum1D(Spectrum1D):
    """
    Subclass of Spectrum1D intended for use in radio astronomy

    Key features (above the Spectrum1D base class): ::
        * Arithmetic operations allow a 'threshold' to be specified, below
          which spectroscopic axes (the "dispersion axis") will be treated as
          identical.  This is implemented because radio telescopes often record
          ~cm/s shifts between scans (because an observatory's relative
          velocity changes with the Earth's rotation)

    """

    def _operation_wrapper(operation):
        """
        Perform an operation (addition, subtraction, mutiplication, division, etc.)
        after checking for shape matching
        """

        def ofunc(self, other): 
            if np.isscalar(other):
                newspec = self.copy()
                newspec.data = operation(newspec.data, other) 
                return newspec
            else: # purely for readability

                if self._arithmetic_threshold == 'exact':
                    dispersioncheck = all(self.dispersion == other.dispersion)
                else:
                    if self._arithmetic_threshold_units is None:
                        # not sure this should ever be allowed
                        dispersioncheck = all((self.dispersion-other.dispersion) < self._arithmetic_threshold)
                    else:
                        dispersioncheck = all((self.dispersion.as_unit(self._arithmetic_threshold_units)-other.dispersion.as_unit(self._arithmetic_threshold_units)) < self._arithmetic_threshold)

                if self.shape == other.shape and dispersioncheck:
                    newspec = self.copy()
                    newspec.data = operation(newspec.data, other.data)
                    return newspec
                elif self.shape != other.shape:
                    raise ValueError("Shape mismatch in data")
                elif not dispersioncheck:
                    raise ValueError("X-axes do not match.")

        return ofunc

    @property
    def _arithmetic_threshold(self):
        return self._arithmetic_threshold_value

    @_arithmetic_threshold.setter
    def _arithmetic_threshold(self, value, units=None):
        self._arithmetic_threshold_value = value
        if units is None:
            self._arithmetic_threshold_units = self.dispersion.units
        else:
            self._arithmetic_threshold_units = units

    _arithmetic_threshold_value = 'exact'
    _arithmetic_threshold_units = None

    __add__ = _operation_wrapper(np.add)
    __radd__ = _operation_wrapper(np.add)
    __sub__ = _operation_wrapper(np.subtract)
    __mul__ = _operation_wrapper(np.multiply)
    __div__ = _operation_wrapper(np.divide)


    # Below are ideas for things to be incorporated into the base Spectrum1D class
    def _shape(self):
        """
        Return the data shape
        """
        return self.data.shape

    shape = property(_shape)

    def __len__(self):
        return len(self.data)

    def __repr__(self):
        # This is a touch verbose, but really nice for interactive use...
        if hasattr(self,'specname'):
            name = " named %s" % self.specname
        else:
            name = ""
        return r'<Spectrum object%s over spectral range %6.5g : %6.5g %s and flux range = [%2.1f, %2.1f] %s at %s>' % \
                (name, self.dispersion.min(), self.dispersion.max(), self.dispersion.units,
                        self.data.min(), self.data.max(), self.units,
                        str(hex(self.__hash__())))
