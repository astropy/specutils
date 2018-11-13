import os

import astropy.units as u
from astropy.nddata import StdDevUncertainty
from astropy.table import Table

from ... import Spectrum1D
from ..registers import data_loader
from ..generic_spectrum_from_table import generic_spectrum_from_table

__all__ = ['ascii_identify', 'ascii_loader', 'ipac_identify', 'ipac_loader']


def ascii_identify(origin, *args, **kwargs):
    """Check if it's an ASCII file."""
    name = os.path.basename(args[0])

    if name.lower().split('.')[-1] in ['txt', 'ascii']:
       return True

    return False


@data_loader(label="ASCII", identifier=ascii_identify, extensions=['txt', 'ascii'])
def ascii_loader(file_name, column_mapping=None, **kwargs):
    """
    Load spectrum from ASCII file.

    Parameters
    ----------
    file_name: str
        The path to the ASCII file
    column_mapping : dict
        A dictionary describing the relation between the ASCII file columns
        and the arguments of the `Spectrum1D` class, along with unit
        information. The dictionary keys should be the ASCII file column names
        while the values should be a two-tuple where the first element is the
        associated `Spectrum1D` keyword argument, and the second element is the
        unit for the ASCII file column::

            column_mapping = {'FLUX': ('flux': 'Jy')}

    Returns
    -------
    data: Spectrum1D
        The spectrum that is represented by the data in this table.
    """
    tab = Table.read(file_name, format='ascii')

    # If not column mapping is given, attempt to parse the ascii files using
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


def ipac_identify(*args, **kwargs):
    """Check if it's an IPAC-style ASCII file."""
    name = os.path.basename(args[0])

    if name.lower().split('.')[-1] in ['txt', 'dat']:
        return True

    return False


@data_loader(label="IPAC", identifier=ipac_identify, extensions=['txt', 'dat'])
def ipac_loader(file_name, column_mapping=None, **kwargs):
    """
    Load spectrum from IPAC-style ASCII file

    Parameters
    ----------
    file_name: str
        The path to the IPAC-style ASCII file.
    column_mapping : dict
        A dictionary describing the relation between the IPAC-style ASCII
        file columns and the arguments of the `Spectrum1D` class, along with
        unit information. The dictionary keys should be the IPAC-style ASCII
        file column names while the values should be a two-tuple where the
        first element is the associated `Spectrum1D` keyword argument, and the
        second element is the unit for the IPAC-style ASCII file column::

            column_mapping = {'FLUX': ('flux': 'Jy')}

    Returns
    -------
    data: Spectrum1D
        The spectrum that is represented by the data in this table.
    """
    tab = Table.read(file_name, format='ascii.ipac')

    # If not column mapping is given, attempt to parse the ascii files using
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
