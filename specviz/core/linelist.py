from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import os

import numpy as np

from astropy.table import Table, vstack

FORMAT = 'line_list'
COLUMN_NAME = 'name'
COLUMN_START = 'start'
COLUMN_END = 'end'
WAVELENGTH_COLUMN = 'wavelength'
ID_COLUMN = 'id'
UNITS_COLUMN = 'units'


# Inheriting from QTable somehow makes this class incompatible
# with the registry machinery in astropy.

class LineList(Table):

    @classmethod
    def ingest(cls, range):
        # Lets skip the file dialog business for now. This is just a
        # proof-of-concept code. Later we will add more fanciness to it.
        #
        # Use these two tables for now. In the future, the filter strings
        # should somehow be handled by the file dialog itself.
        fnames = ['Common_stellar.txt', 'Common_nebular.txt']

        path = os.path.dirname(os.path.abspath(__file__))
        dir_path = path + '/../data/linelists/'
        linelists = []

        for fname in fnames:
            path = dir_path + fname
            filter = fname.split('.')[0] + ' (*.txt *.dat)'

            linelist = LineList.read(path, filter)
            linelist = linelist.extract_range(range)

            linelists.append(linelist)

        return LineList.merge(linelists)

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

    def extract_range(self, wrange):
        ''' Builds a LineList instance out of self, with
            the subset of lines that fall within the
            wavelength range defined by 'wmin' and 'wmax'

        :param wrange: tuple of 2 floats
            minimum and maximum wavelength of the wavelength range
        :return: LineList
            line list with subset of lines
        '''
        wavelengths = self[WAVELENGTH_COLUMN].quantity

        wmin = wrange[0]
        wmax = wrange[1]

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
