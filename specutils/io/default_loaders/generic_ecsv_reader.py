import os

from astropy.nddata import StdDevUncertainty
from astropy.table import Table
import astropy.units as u

from ...spectra import Spectrum
from ..registers import data_loader, custom_writer
from ..parsing_utils import (generic_spectrum_from_table,
                             spectrum_from_column_mapping)


def identify_ecsv(origin, *args, **kwargs):
    """Check if it's an ECSV file."""
    return (isinstance(args[0], (str, os.PathLike)) and
            os.path.splitext(os.fspath(args[0]).lower())[1] == '.ecsv')


@data_loader("ECSV", identifier=identify_ecsv, dtype=Spectrum, priority=-10)
def generic_ecsv(file_name, column_mapping=None, **kwargs):
    """
    Read a spectrum from an ECSV file, using generic_spectrum_from_table_loader()
    to try to figure out which column is which.
    The ECSV columns must have units, as `generic_spectrum_from_table_loader`
    depends on this to determine the meaning of the columns.  For manual
    control over the column to spectrum mapping, use the ASCII loader.

    Parameters
    ----------
    file_name: str
        The path to the ECSV file.
    column_mapping : dict
        A dictionary describing the relation between the ECSV file columns
        and the arguments of the `Spectrum` class, along with unit
        information. The dictionary keys should be the ECSV file column names
        while the values should be a two-tuple where the first element is the
        associated `Spectrum` keyword argument, and the second element is the
        unit for the ECSV file column::

            column_mapping = {'FLUX': ('flux', 'Jy')}

    Returns
    -------
    data: Spectrum
        The spectrum that is represented by the data in this table.
    """
    table = Table.read(file_name, format='ascii.ecsv')

    if column_mapping is None:
        return generic_spectrum_from_table(table)

    return spectrum_from_column_mapping(table, column_mapping)


@custom_writer("ECSV")
def generic_ecsv_writer(spectrum, file_name, **kwargs):
    """
    Write a spectrum to an ECSV table.

    Parameters
    ----------
    spectrum : :class:`Spectrum`
        The spectrum to write.
    file_name : str, file-like or `pathlib.Path`
        File to write to.
    **kwargs
        Additional optional keywords passed to :meth:`astropy.table.Table.write`.
    """
    columns = [spectrum.spectral_axis, spectrum.flux]
    colnames = ["spectral_axis", "flux"]

    # The generic ECSV reader interprets a second flux-like column as a
    # standard-deviation uncertainty, matching the tabular FITS convention.
    if spectrum.uncertainty is not None:
        uncertainty = (
            spectrum.uncertainty
            .represent_as(StdDevUncertainty)
            .quantity
            .to(spectrum.flux.unit,
                equivalencies=u.spectral_density(spectrum.spectral_axis))
        )
        columns.append(uncertainty)
        colnames.append("uncertainty")

    if spectrum.mask is not None:
        columns.append(spectrum.mask)
        colnames.append("mask")

    # Astropy tables store the spectral coordinate along the row axis. Match
    # the tabular FITS writer by moving the spectral dimension to the front for
    # multidimensional flux, uncertainty, and mask arrays.
    for index in range(1, len(columns)):
        if columns[index].ndim > 1:
            columns[index] = columns[index].T

    # ECSV metadata must be YAML-serializable. Converting a FITS Header to a
    # plain dict preserves its values while matching the reader's convention of
    # storing table metadata in Spectrum.meta["header"].
    meta = dict(spectrum.meta.get("header", {}))
    table = Table(columns, names=colnames, meta=meta)
    table.write(file_name, format="ascii.ecsv", **kwargs)
