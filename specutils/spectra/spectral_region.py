import numpy as np
import astropy.units as u


class SpectralRegion:
    """
    A `SpectralRegion` is a container class enables some simplicty
    in defining and passing a region (interval) for a spectrum.

    In the future, there might be more functionality added in here and there
    is some discussion that this might/could move to
    `Astropy Regions <http://astropy-regions.readthedocs.io/en/latest/>`_.
    """

    @staticmethod
    def from_center(center=None, width=None):
        """
        SpectralRegion class method that enables the definition of a `SpectralRegion`
        from the center and width rather than lower and upper bounds.

        Parameters
        ----------

        center: Scalar `~astropy.units.Quantity` with pixel or ``spectral_axis`` units
           The center of the spectral region.

        Width: Scalar `~astropy.units.Quantity` with pixel or ``spectral_axis`` units
           The width of the spectral region.
        """

        return SpectralRegion(center - width, center + width)

    def __init__(self, lower, upper):
        """
        Lower and upper values for the interval.

        Parameters
        ----------

        lower: Scalar `~astropy.units.Quantity` with pixel or ``spectral_axis`` units
           The lower bound of the region.

        upper: Scalar `~astropy.units.Quantity` with pixel or ``spectral_axis`` units
           The upper bound of the region.

        """

        # Create instance variables
        self._lower = None
        self._upper = None

        # Set the values (using the setters for doing the proper checking)
        self.lower = lower
        self.upper = upper

    @property
    def lower(self):
        return self._lower

    @lower.setter
    def lower(self, value):
        """
        Lower bound for the interval.

        Parameters
        ----------

        value: Scalar `~astropy.units.Quantity` with pixel or ``spectral_axis`` units
           The lower bound of the region.
        """

        if not (value.unit.is_equivalent(u.pixel) or
           value.unit.is_equivalent(u.angstrom, equivalencies=u.spectral())):
            raise u.UnitError('Lower bound of region is not a spectral region unit')

        self._lower = value

    @property
    def upper(self):
        return self._upper

    @upper.setter
    def upper(self, value):
        """
        Upper bound for the interval.

        Parameters
        ----------

        value: Scalar `~astropy.units.Quantity` with pixel or ``spectral_axis`` units
           The upper bound of the region.
        """

        if not (value.unit.is_equivalent(u.pixel) or
           value.unit.is_equivalent(u.angstrom, equivalencies=u.spectral())):
            raise u.UnitError('Upper bound of region is not a spectral region unit')

        self._upper = value

    def to_pixel(self, spectrum):
        """
        Calculate and return the left and right indices defined
        by the lower and upper bounds and based on the input 
        `~specutils.spectra.spectrum1d.Spectrum1D`. The left and right indices will
        be returned.

        Parameters
        ----------
        spectrum: `~specutils.spectra.spectrum1d.Spectrum1D`
            The spectrum object from which the region will be extracted.

        Return
        ------
        left_index, right_index: int, int
            Left and right indices defined by the lower and upper bounds.

        """

        left_index = int(np.ceil(spectrum.wcs.world_to_pixel(np.array(self._lower))))
        right_index = int(np.floor(spectrum.wcs.world_to_pixel(np.array(self._upper))))

        return left_index, right_index

    def extract(self, spectrum):
        """
        Extract a region from the input `~specutils.spectra.spectrum1d.Spectrum1D`
        defined by the `lower` and `upper` bounds defined by this SpectralRegion
        instance.  The extracted region will be returned.

        Parameters
        ----------
        spectrum: `~specutils.spectra.spectrum1d.Spectrum1D`
            The spectrum object from which the region will be extracted.

        Return
        ------
        spectrum: `~specutils.spectra.spectrum1d.Spectrum1D`
            Excised spectrum.

        Notes
        -----
        The region extracted is a discrete subset of the input spectrum. No interpolation is done
        on the left and right side of the spectrum.

        """

        left_index, right_index = self.to_pixel(spectrum)

        if left_index >= right_index:
            raise ValueError('Lower region, {}, appears to be greater than the upper region, {}.'.format(self._lower, self._upper))

        return spectrum[left_index:right_index]
