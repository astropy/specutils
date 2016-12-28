"""
Emission/Absorption Line list utilities
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import os
import glob

import numpy as np

from astropy.table import Table, vstack

__all__ = [
    'LineList',
    'ingest',
]

FORMAT = 'line_list'
COLUMN_NAME = 'name'
COLUMN_START = 'start'
COLUMN_END = 'end'
WAVELENGTH_COLUMN = 'Wavelength'
ID_COLUMN = 'Line ID'
UNITS_COLUMN = 'units'

def ingest(range):
    """
    Returns a list with LineList instances.

    Each original list is stripped out of lines that lie outside the
    wavelength range.

    Parameters
    ----------
    range:
        The wavelength range of interest.

    Returns
    -------
    [LineList, ...]
        The list of linelists found.

    Notes
    -----
    Lets skip the file dialog business. For now, we look
    for line lists and their accompanying YAML files in
    one single place. We also restrict our search for
    ascii line lists whose file names end in .txt
    """
    linelist_path = os.path.dirname(os.path.abspath(__file__))
    dir_path = linelist_path + '/../data/linelists/'
    yaml_paths = glob.glob(dir_path + '*.yaml')
    linelists = []

    for yaml_path in yaml_paths:
        # this should get improved as when we decide how to
        # implement support for user-supplied line lists,
        # as well as support for other formats besides ascii.
        linelist_path = yaml_path.replace('.yaml', '.txt')
        filter = linelist_path.split(os.sep)[-1].split('.')[0] + ' (*.txt *.dat)'

        linelist = LineList.read(linelist_path, format="ascii.tab")
        linelist = linelist.extract_range(range)

        linelists.append(linelist)

    return linelists


# Inheriting from QTable somehow makes this class incompatible
# with the registry machinery in astropy.

class LineList(Table):
    """
    A list of emission/absorption lines

    Parameters
    ----------
    table: `~astropy.table.Table`
        If specified, a table to initialize from.

    name: str
        The name of the list.

    masked: bool
        If true, a masked table is used.
    """

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
        """
        Executes a 'vstack' of all input lists, and
        then sorts the result by the wavelength column.

        Parameters
        ----------
        lists: [LineList, ...]
            list of LineList instances

        Returns
        -------
        LineList
            merged line list
        """
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
        """
        Builds a LineList instance out of self, with
        the subset of lines that fall within the
        wavelength range defined by 'wmin' and 'wmax'

        Parameters
        ----------
        wrange: (float, float)
            minimum and maximum wavelength of the wavelength range

        Returns
        -------
        LineList
            line list with subset of lines
        """
        wavelengths = self[WAVELENGTH_COLUMN].quantity

        wmin = wrange[0]
        wmax = wrange[1]

        # convert wavelenghts in line list to whatever
        # units the wavelength range is expressed in.
        new_wavelengths = wavelengths.to(wmin.unit)

        # 'indices' points to rows with wavelength values
        # that lie outside the wavelength range.
        indices_to_remove = np.where((new_wavelengths.value < wmin.value) |
                                     (new_wavelengths.value > wmax.value))

        return self._remove_lines(indices_to_remove)

    def extract_rows(self, indices):
        """
        Builds a LineList instance out of self, with
        the subset of lines pointed by 'indices'

        Parameters
        ----------
        indices: [QModelIndex, ...]
            List of QModelIndex instances to extract from.

        Returns
        -------
        LineList
            line list with subset of lines
        """
        row_indices = []
        for index in indices:
            row_indices.append(index.row())

        line_indices = []
        for index in range(len(self.columns[0])):
            line_indices.append(index)

        indices_to_remove = list(
            filter(lambda x: x not in row_indices, line_indices)
        )

        return self._remove_lines(indices_to_remove)

    def _remove_lines(self, indices_to_remove):
        """
        Makes a copy of self and removes
        unwanted lines from the copy.

        Parameters
        ----------
        indices_to_remove: [int, ...]
            List of row numbers to remove

        Returns
        -------
        LineList:
            A new copy of the `LineList` with the rows removed.
        """
        table = Table(self)

        table.remove_rows(indices_to_remove)

        result = LineList(table, self.name)

        return result
