from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from astropy.table import Table

COLUMN_NAME = 'name'
COLUMN_START = 'start'
COLUMN_END = 'end'
WAVELENGTH_COLUMN = 'wavelength'
ID_COLUMN = 'ids'
UNITS_COLUMN = 'units'

class LineList(Table):

    @classmethod
    def read(cls, *args, **kwargs):
        from ..interfaces.registries import io_registry
        return io_registry.read(cls, *args, **kwargs)

    def extract_range(self, wmin, wmax):
        wavelengths = self[WAVELENGTH_COLUMN].quantity
        new_wavelengths = wavelengths.to(wmin.unit)

        indices = np.where(new_wavelengths.value.any() >= wmin.value and
                           new_wavelengths.value.any() <= wmax.value)

        print("@@@@@@  file linelist.py; line 26 - ",  wmin, wmax, new_wavelengths.unit, indices)
        print("@@@@@@  file linelist.py; line 29 - ",  new_wavelengths)
