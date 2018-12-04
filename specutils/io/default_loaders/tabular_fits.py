import logging
import os

import numpy as np

from astropy.io import fits
from astropy.nddata import StdDevUncertainty
from astropy.table import Table
import astropy.units as u
from astropy.wcs import WCS

from ...spectra import Spectrum1D
from ..registers import data_loader, custom_writer
from ..parsing_utils import (generic_spectrum_from_table,
                             spectrum_from_column_mapping)

__all__ = ['tabular_fits_loader', 'tabular_fits_writer']


def identify_tabular_fits(origin, *args, **kwargs):
    # check if file can be opened with this reader
    # args[0] = filename
    # fits.open(args[0]) = hdulist
    return (isinstance(args[0], str) and
            # check if file is .fits
            os.path.splitext(args[0].lower())[1] == '.fits' and
            # check hdulist has more than one extension
            len(fits.open(args[0])) > 1 and
            # check if fits has BinTable extension
            isinstance(fits.open(args[0])[1], fits.BinTableHDU)
            )


@data_loader("tabular-fits", identifier=identify_tabular_fits,
             dtype=Spectrum1D, extensions=['fits'])
def tabular_fits_loader(file_name, column_mapping=None, **kwargs):
    """
    Load spectrum from a FITS file.

    Parameters
    ----------
    file_name: str
        The path to the FITS file
    column_mapping : dict
        A dictionary describing the relation between the FITS file columns
        and the arguments of the `Spectrum1D` class, along with unit
        information. The dictionary keys should be the FITS file column names
        while the values should be a two-tuple where the first element is the
        associated `Spectrum1D` keyword argument, and the second element is the
        unit for the ASCII file column::

            column_mapping = {'FLUX': ('flux', 'Jy')}

    Returns
    -------
    data: Spectrum1D
        The spectrum that is represented by the data in this table.
    """
    # Parse the wcs information. The wcs will be passed to the column finding
    # routines to search for spectral axis information in the file.
    with fits.open(file_name) as hdulist:
        wcs = WCS(hdulist[0].header)

    tab = Table.read(file_name, format='fits')

    # If no column mapping is given, attempt to parse the file using
    # unit information
    if column_mapping is None:
        return generic_spectrum_from_table(tab, wcs=wcs, **kwargs)

    return spectrum_from_column_mapping(tab, column_mapping, wcs=wcs)


@custom_writer("tabular-fits")
def tabular_fits_writer(spectrum, file_name, **kwargs):
    flux = spectrum.flux.value
    disp = spectrum.dispersion.value
    meta = spectrum.meta

    tab = Table([disp, flux], names=("dispersion", "flux"), meta=meta)

    tab.write(file_name, format="fits")
