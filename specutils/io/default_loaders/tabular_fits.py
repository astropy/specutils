import logging
import os

import six
from astropy.io import fits
from astropy.nddata import StdDevUncertainty
from astropy.table import Table
import astropy.units as u
import numpy as np

from ...spectra import Spectrum1D
from ..registers import data_loader, custom_writer
from ..generic_spectrum_from_table import generic_spectrum_from_table

__all__ = ['tabular_fits_loader', 'tabular_fits_writer']


def identify_tabular_fits(origin, *args, **kwargs):
    # check if file can be opened with this reader
    # args[0] = filename
    # fits.open(args[0]) = hdulist
    return (isinstance(args[0], six.string_types) and
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

            column_mapping = {'FLUX': ('flux': 'Jy')}

    Returns
    -------
    data: Spectrum1D
        The spectrum that is represented by the data in this table.
    """
    tab = Table.read(file_name, format='fits')

    # If not column mapping is given, attempt to parse the file using
    # unit information
    if column_mapping is None:
        return generic_spectrum_from_table(tab, **kwargs)

    spec_kwargs = {}

    # Associate columns of the file with the appropriate spectrum1d arguments
    for col_name, (kwarg_name, unit) in column_mapping.items():
        if unit is None:
            kwarg_val = tab[column_mapping]
        else:
            kwarg_val = u.Quantity(tab[column_mapping], unit)

        spec_kwargs.setdefault(kwarg_name, kwarg_val)

    # Ensure that the uncertainties are a subclass of NDUncertainty
    if spec_kwargs.get('uncertainty') is not None:
        spec_kwargs['uncertainty'] = StdDevUncertainty(spec_kwargs.get('uncertainty'))

    return Spectrum1D(**spec_kwargs, meta=tab.meta)


@custom_writer("tabular-fits")
def tabular_fits_writer(spectrum, file_name, **kwargs):
    flux = spectrum.flux.value
    disp = spectrum.dispersion.value
    meta = spectrum.meta

    tab = Table([disp, flux], names=("dispersion", "flux"), meta=meta)

    tab.write(file_name, format="fits")