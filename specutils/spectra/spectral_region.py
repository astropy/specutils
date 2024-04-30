import itertools
import sys

from astropy.table import QTable
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
        self._valid()

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

        if center.unit.physical_type not in ('length', 'unknown'):
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

    @classmethod
    def from_qtable(cls, table):
        """
        Generate a ``SpectralRegion`` instance from an `~astropy.table.QTable`
        object has has ``lower_bound`` and ``upper_bound`` columns

        Parameters
        ----------
        table : `~astropy.table.QTable`
            An `~astropy.table.QTable` object with ``lower_bound`` and ``upper_bound`` columns.

        Returns
        -------
        `~specutils.SpectralRegion`
            The spectral region based on the table of bounds.
        """
        subregions = []
        for row in table:
            subregions.append([row["lower_bound"], row["upper_bound"]])

        return cls(subregions)

    @classmethod
    def read(cls, filename):
        """
        Create a ``SpectralRegion`` from an ecsv file output by the
        ``SpectralRegion.write`` method.

        Parameters
        ----------
        filename : str
            The name of the ecsv file on disk to be read in as a ``SpectralRegion``.
        """
        table = QTable.read(filename)
        return cls.from_qtable(table)

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

        bound_unit = self._subregions[0][0].unit
        for x in self._subregions:
            if x[0].unit != bound_unit or x[1].unit != bound_unit:
                raise ValueError("All SpectralRegion bounds must have the same unit.")
            if x[0] == x[1]:
                raise ValueError("Upper and lower bound must be different values.")

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

    @property
    def _table(self):
        # Return the information defining the spectral region as an astropy QTable
        lower_bounds = []
        upper_bounds = []
        for subregion in self.subregions:
            lower_bounds.append(subregion[0])
            upper_bounds.append(subregion[1])

        return QTable([lower_bounds, upper_bounds], names=("lower_bound", "upper_bound"))

    def as_table(self):
        """Returns an `~astropy.table.QTable` with the upper and lower bound
        of each subregion in the ``SpectralRegion``."""
        return self._table

    def write(self, filename="spectral_region.ecsv", overwrite=True):
        """
        Write the SpectralRegion to an ecsv file using `~astropy.table.QTable`.
        Overwrites by default.
        """
        if not filename.endswith("ecsv"):
            raise ValueError("SpectralRegions can only be written out to ecsv files.")

        self._table.write(filename, overwrite=overwrite)
