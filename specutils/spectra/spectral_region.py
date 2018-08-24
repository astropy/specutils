import itertools
import sys

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

        center: Scalar `~astropy.units.Quantity` with pixel or any valid ``spectral_axis`` unit
           The center of the spectral region.

        Width: Scalar `~astropy.units.Quantity` with pixel or any valid ``spectral_axis`` unit
           The width of the spectral region.
        """

        return SpectralRegion(center - width, center + width)

    def __init__(self, *args):
        """
        Lower and upper values for the interval.

        Parameters
        ----------

        lower: Scalar `~astropy.units.Quantity` with pixel or any valid ``spectral_axis`` unit
           The lower bound of the region.

        upper: Scalar `~astropy.units.Quantity` with pixel or any valid ``spectral_axis`` unit
           The upper bound of the region.

        """

        # Create instance variables
        self._subregions = None

        #
        # Set the values (using the setters for doing the proper checking)
        #

        if self._is_2_element(args):
            self._subregions = [tuple(args)]
        elif isinstance(args, (list, tuple)) and all([self._is_2_element(x) for x in args[0]]):
            self._subregions = [tuple(x) for x in args[0]]
        else:
            raise ValueError('SpectralRegion input must be a 2-tuple or a list of 2-tuples.')

        #
        #  Check validity of the input sub regions.
        #
        if not self._valid():
            raise ValueError('SpectralRegion 2-tuple lower extent must be less than upper extent.')

    def __str__(self):
        return 'SpectralRegion: {}'.format(
                ', '.join(['{} - {}'.format(x[0], x[1]) for x in self._subregions]))

    def __add__(self, other):
        """
        Ability to add two SpectralRegion classes together.
        """
        return SpectralRegion(self._subregions + other._subregions)

    def __iadd__(self, other):
        """
        Ability to add one SpectralRegion to another using +=.
        """
        self._subregions += other._subregions
        return self

    def __len__(self):
        """
        Number of spectral regions.
        """
        return len(self._subregions)

    def __getslice__(self, item):
        """
        Enable slicing of the SpectralRegion list.
        """
        return SpectralRegion(self._subregions[item])

    def __getitem__(self, item):
        """
        Enable slicing or extracting the SpectralRegion.
        """
        if isinstance(item, slice):
            return self.__getslice__(item)
        else:
            return SpectralRegion([self._subregions[item]])

    def __delitem__(self, item):
        """
        Delete a specific item from the list.
        """
        del self._subregions[item]

    def _valid(self):
        return all(x[0] < x[1] for x in self._subregions)

    def _is_2_element(self, value):
        """
        Helper function to check a variable to see if it
        is a 2-tuple.
        """
        return len(value) == 2 and \
               isinstance(value[0], u.Quantity) and \
               isinstance(value[1], u.Quantity)

    def _reorder(self):
        """
        Re-order the  list based on lower bounds.

        This could be important for when a spectrum is extracted.
        """
        self._subregions.sort(key=lambda k: k[0])

    @property
    def bounds(self):
        """
        Compute the lower and upper extent of the SpectralRegion.
        """
        return (self.lower, self.upper)

    @property
    def lower(self):
        """
        The most minimum value of the sub-regions.
        """
        return min(x[0] for x in self._subregions)

    @property
    def upper(self):
        """
        The most maximum value of the sub-regions.
        """
        return max(x[0] for x in self._subregions)

    def invert_from_spectrum(self, spectrum):
        """
        Invert a SpectralRegion based on the extent of the
        input spectrum.  

        See notes in SpectralRegion.invert() method.
        """
        return self.invert(spectrum.spectral_axis[0], spectrum.spetral_axis[-1])

    def _in_range(self, value, lower, upper):
        return (value >= lower) and (value <= upper)

    def invert(self, lower_bound, upper_bound):
        """
        Given a set of sub-regions this SpectralRegion defines,
        create a new SpectralRegion such that the sub-regions
        are defined in the new one as regions *not* in this
        SpectralRegion.

        For example, if this SpectralRegion is defined as:
           sr = SpectralRegion([(0.45*u.um, 0.6*u.um), (0.8*u.um, 0.9*u.um)])

        and we call `sr_invert = sr.invert(0.3*u.um, 1.0*u.um)` then
        `sr_invert` will be
           SpectralRegion([(0.3*u.um, 0.45*u.um), (0.6*u.um, 0.8*u.um), (0.9*u.um, 1*u.um)])

        This could be useful if a SpectralRegion has sub-regions defined
        for peaks in a spectrum and then one wants to create a SpectralRegion
        defined as all the *non*-peaks, then one could use this function.
        """

        #
        # Create 'rs' region list with left and right extra ranges.
        #
        min_num = -sys.maxsize-1
        max_num = sys.maxsize
        rs = self._subregions + [(min_num*u.um, lower_bound), (upper_bound, max_num*u.um)]

        #
        # Sort the region list based on lower bound.
        #

        sorted_regions = sorted(rs, key=lambda k: k[0])

        #
        # Create new region list that has overlapping regions merged
        #

        merged = []
        for higher in sorted_regions:
            if not merged:
                merged.append(higher)
            else:
                lower = merged[-1]
                # test for intersection between lower and higher:
                # we know via sorting that lower[0] <= higher[0]
                if higher[0] <= lower[1]:
                    upper_bound = max(lower[1], higher[1])
                    merged[-1] = (lower[0], upper_bound)  # replace by merged interval
                else:
                    merged.append(higher)

        #
        # Create new list and drop first and last (the maxsize ones).
        # We go from -inf, upper1, lower2, upper2....
        # and remap to     lower1, upper1, lower2, ...
        #

        newlist = list(itertools.chain.from_iterable(merged))
        newlist = newlist[1:-1]

        #
        # Now create new Spectrum region
        #

        return SpectralRegion([(x, y) for x, y in zip(newlist[0::2], newlist[1::2])])

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

        left_index = int(np.ceil(spectrum.wcs.world_to_pixel(self.lower)))
        right_index = int(np.floor(spectrum.wcs.world_to_pixel(self.upper)))

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
