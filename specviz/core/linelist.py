from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import os
import glob

import numpy as np

from astropy.table import Table, vstack

FORMAT = 'line_list'
COLUMN_NAME = 'name'
COLUMN_START = 'start'
COLUMN_END = 'end'
WAVELENGTH_COLUMN = 'Wavelength'
ID_COLUMN = 'Line ID'
UNITS_COLUMN = 'units'


# Returns a list with LineList instances. Each original list is
# stripped out of lines that lie outside the wavelength range.

def ingest(range):

    # Lets skip the file dialog business. For now, we look
    # for line lists and their accompanying YAML files in
    # one single place. We also restrict our search for
    # ascii line lists whose file names end in .txt

    path = os.path.dirname(os.path.abspath(__file__))
    dir_path = path + '/../data/linelists/'
    yaml_paths = glob.glob(dir_path + '*.yaml')
    linelists = []

    for yaml_path in yaml_paths:
        # this should get improved as when we decide how to
        # implement support for user-supplied line lists,
        # as well as support for other formats besides ascii.
        path = yaml_path.replace('.yaml', '.txt')
        filter = path.split(os.sep)[-1].split('.')[0] + ' (*.txt *.dat)'

        linelist = LineList.read(path, filter)
        linelist = linelist.extract_range(range)

        linelists.append(linelist)

    return linelists


# Inheriting from QTable somehow makes this class incompatible
# with the registry machinery in astropy.

class LineList(Table):

    def __init__(self, table=None, name=None, masked=None):
        Table.__init__(self, data=table, masked=masked)

        self.name = name

        # We have to carry internally a raw reference to the
        # table data so as to be able to use vstack() to perform
        # merging. This shouldn't be a problem as long as the
        # LineList instance is regarded as immutable. Which it
        # should be anyways.

        self._table = table

    @classmethod
    def merge(cls, lists):
        ''' Executes a 'vstack' of all input lists, and
            then sorts the result by the wavelength column.

        :param lists: list
            list of LineList instances
        :return: LineList
            merged line list
        '''
        # Note that vstack operates on Table instances but
        # not on LineList instances. So we first extract the
        # raw Table instances.
        tables = []
        for linelist in lists:
            tables.append(linelist._table)

        merged_table = vstack(tables)

        merged_table.sort(WAVELENGTH_COLUMN)

        return cls(merged_table, "Merged")

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

        result = LineList(result, self.name)

        return result
