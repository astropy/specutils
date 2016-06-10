from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from astropy.table import Table, vstack

FORMAT = 'line_list'
COLUMN_NAME = 'name'
COLUMN_START = 'start'
COLUMN_END = 'end'
WAVELENGTH_COLUMN = 'wavelength'
ID_COLUMN = 'ids'
UNITS_COLUMN = 'units'

# Inheriting from QTable somehow makes this class incompatible
# with the registry machinery in astropy.

class LineList(Table):

    @classmethod
    def merge(cls, lists):
        ''' Executes a 'vstack' of all input lists, and
            then sorts the result by the wavelength column.

        :param lists: list
            list of LineList instances
        :return: LineList
            merged line list
        '''
        result = vstack(lists)

        result.sort(WAVELENGTH_COLUMN)

        return result

    def extract_range(self, wmin, wmax):
        ''' Builds a LineList instance out of self, with
            the subset of lines that fall within the
            wavelength range defined by 'wmin' and 'wmax'

        :param wmin: float
            minimum wavelength of the wavelength range
        :param wmax: float
            maximum wavelength of the wavelength range
        :return: LineList
            line list with subset of lines
        '''
        wavelengths = self[WAVELENGTH_COLUMN].quantity

        # convert wavelenghts in line list to whatever
        # units the wavelength range is expressed in.
        new_wavelengths = wavelengths.to(wmin.unit)

        # 'indices' points to rows with wavelength values
        # that lie outside the wavelength range.
        indices = np.where((new_wavelengths.value < wmin.value) |
                           (new_wavelengths.value > wmax.value))

        # make copy of self and remove unwanted lines from the copy.
        result = Table(self)
        result.remove_rows(indices)

        return result
