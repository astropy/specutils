import itertools
import sys

import numpy as np
import astropy.units as u


class SpectralRegion:
    """
    A `SpectralRegion` is a container class for regions (intervals) along a
    spectral coordinate.

    This class can either represent a single contiguous region or a set of
    regions related to each other in some way (For example, a pair of
    continuum windows around a line or a doublet of lines).

    Parameters
    ----------
    *args : variable
        Either a single parameter or two parameters can be given:

        * 1 argument ``regioniter``: An iterable of length-2
          `~astropy.units.Quantity` objects (or a single n x 2
          `~astropy.units.Quantity` object), where the length-2 dimension is
          ``lower``, ``upper``.
        * ``lower``, ``upper``: Each should be an `~astropy.units.Quantity`
          object.

    Notes
    -----
    The subregions will be ordered based on the lower bound of each subregion.

    """
    def __init__(self, *args):
        # Create instance variables
        self._subregions = None

        # Set the values (using the setters for doing the proper checking)
        if len(args) == 1:
            if all([self._is_2_element(x) for x in args[0]]):
                self._subregions = [tuple(x) for x in args[0]]
            else:
                raise ValueError("SpectralRegion 1-argument input must be "
                                 "a list of length-two Quantity objects.")
        elif len(args) == 2:
            if self._is_2_element(args):
                self._subregions = [tuple(args)]
            else:
                raise ValueError("SpectralRegion 2-argument inputs must be "
                                 "Quantity objects.")
        else:
            raise TypeError(f'SpectralRegion initializer takes 1 or 2'
                            f'positional arguments but {len(args)} were given')

        #  Check validity of the input sub regions.
        if not self._valid():
            raise ValueError("SpectralRegion 2-tuple lower extent must be "
                             "less than upper extent.")

        # The sub-regions are to always be ordered based on the lower bound.
        self._reorder()

    @classmethod
    def from_center(cls, center=None, width=None):
        """
        SpectralRegion class method that enables the definition of a
        `SpectralRegion` from the center and width rather than lower and
        upper bounds.

        Parameters
        ----------
        center : Scalar `~astropy.units.Quantity` with pixel or any valid ``spectral_axis`` unit
           The center of the spectral region.
        width : Scalar `~astropy.units.Quantity` with pixel or any valid ``spectral_axis`` unit
           The full width of the spectral region (upper bound - lower bound).
        """

        if width.value <= 0:
            raise ValueError("SpectralRegion width must be positive.")

        if center.unit.physical_type != 'length':
            return cls(center + width/2, center - width/2)

        return cls(center - width/2, center + width/2)

    @classmethod
    def from_line_list(cls, table, width=1):
        """
        Generate a ``SpectralRegion`` instance from the `~astropy.table.QTable`
        object returned from `~specutils.fitting.find_lines_derivative` or
        `~specutils.fitting.find_lines_threshold`.

        Parameters
        ----------
        table : `~astropy.table.QTable`
            List of found lines.
        width : float
            The width of the spectral line region. If not unit information is
            provided, it's assumed to be the same units as used in the line
            list table.

        Returns
        -------
        `~specutils.SpectralRegion`
            The spectral region based on the line list.
        """
        width = u.Quantity(width, table['line_center'].unit)

        return cls([(x - width * 0.5, x + width * 0.5)
                    for x in table['line_center']])

    def _info(self):
        """
        Pretty print the sub-regions.
        """

        toreturn = "Spectral Region, {} sub-regions:\n".format(
            len(self._subregions))

        # Setup the subregion text.
        subregion_text = []
        for ii, subregion in enumerate(self._subregions):
            subregion_text.append('  ({}, {})'.format(subregion[0], subregion[1]))

        # Determine the length of the text boxes.
        max_len = max(len(srt) for srt in subregion_text) + 1
        ncols = 70 // max_len

        # Add sub region info to the output text.
        fmt = '{' + ':<{}'.format(max_len) + '}'
        for ii, srt in enumerate(subregion_text):
            toreturn += fmt.format(srt)
            if ii % ncols == (ncols-1):
                toreturn += '\n'

        return toreturn

    def __str__(self):
        return self._info()

    def __repr__(self):
        return self._info()

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
        self._reorder()
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

        # Lower bound < Upper bound for all sub regions in length physical type
        with u.set_enabled_equivalencies(u.spectral()):
            sub_regions = [(x[0].to('m'), x[1].to('m'))
                           if x[0].unit.is_equivalent(u.m) else x
                           for x in self._subregions]

        if any(x[0] >= x[1] for x in sub_regions):
            raise ValueError('Lower bound must be strictly less than the upper bound')

        return True

    @staticmethod
    def _is_2_element(value):
        """
        Helper function to check a variable to see if it
        is a 2-tuple of Quantity objects.
        """
        return len(value) == 2 and \
               isinstance(value[0], u.Quantity) and \
               isinstance(value[1], u.Quantity)

    def _reorder(self):
        """
        Re-order the  list based on lower bounds.
        """
        self._subregions.sort(key=lambda k: k[0])

    @property
    def subregions(self):
        """
        An iterable over ``(lower, upper)`` tuples that are each of the
        sub-regions.
        """
        return self._subregions

    @property
    def bounds(self):
        """
        Compute the lower and upper extent of the SpectralRegion.
        """
        return self.lower, self.upper

    @property
    def lower(self):
        """
        The most minimum value of the sub-regions.

        The sub-regions are ordered based on the lower bound, so the
        lower bound for this instance is the lower bound of the first
        sub-region.
        """
        return self._subregions[0][0]

    @property
    def upper(self):
        """
        The most maximum value of the sub-regions.

        The sub-regions are ordered based on the lower bound, but the
        upper bound might not be the upper bound of the last sub-region
        so we have to look for it.
        """
        return max(x[1] for x in self._subregions)

    def invert_from_spectrum(self, spectrum):
        """
        Invert a SpectralRegion based on the extent of the
        input spectrum.

        See notes in SpectralRegion.invert() method.
        """
        return self.invert(spectrum.spectral_axis[0],
                           spectrum.spectral_axis[-1])

    def _in_range(self, value, lower, upper):
        return (value >= lower) and (value <= upper)

    def invert(self, lower_bound, upper_bound):
        """
        Invert this spectral region.  That is, given a set of sub-regions this
        object defines, create a new `SpectralRegion` such that the sub-regions
        are defined in the new one as regions *not* in this `SpectralRegion`.

        Parameters
        ----------
        lower_bound : `~astropy.units.Quantity`
           The lower bound of the region. Can be scalar with pixel or any
           valid ``spectral_axis`` unit
        upper_bound : `~astropy.units.Quantity`
           The upper bound of the region. Can be scalar with pixel or any
           valid ``spectral_axis`` unit

        Returns
        -------
        spectral_region : `~specutils.SpectralRegion`
           Spectral region of the non-selected regions

        Notes
        -----
        This is applicable if, for example, a `SpectralRegion` has sub-regions
        defined for peaks in a spectrum and then one wants to create a
        `SpectralRegion` defined as all the *non*-peaks, then one could use this
        function.

        As an example, assume this SpectralRegion is defined as
        ``sr = SpectralRegion([(0.45*u.um, 0.6*u.um), (0.8*u.um, 0.9*u.um)])``.
        If we call ``sr_invert = sr.invert(0.3*u.um, 1.0*u.um)`` then
        ``sr_invert`` will be
        ``SpectralRegion([(0.3*u.um, 0.45*u.um), (0.6*u.um, 0.8*u.um), (0.9*u.um, 1*u.um)])``

        """

        #
        # Create 'rs' region list with left and right extra ranges.
        #
        min_num = -sys.maxsize-1
        max_num = sys.maxsize
        rs = self._subregions + [(min_num*u.um, lower_bound),
                                 (upper_bound, max_num*u.um)]

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
