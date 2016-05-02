"""This module contains functions that perform the actual data parsing."""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# STDLIB
import logging
import os

# THIRD-PARTY
import numpy as np
from astropy import units as u
from astropy.io import ascii, fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.nddata import StdDevUncertainty

# LOCAL
from ..core.data import Data
from ..core.linelist import LineList

__all__ = ['fits_reader', 'fits_identify',
           'ascii_reader', 'ascii_identify',
           'linelist_reader', 'linelist_identify']

# Loader automatically falls back to these units for some cases
default_waveunit = u.Unit('Angstrom')
default_fluxunit = u.Unit('erg / (Angstrom cm2 s)')


def fits_reader(filename, filter, **kwargs):
    """This generic function will query the loader factory, which has already
    loaded the YAML configuration files, in an attempt to parse the
    associated FITS file.

    Parameters
    ----------
    filename : str
        Input filename.

    filter : str
        File type for YAML look-up.

    kwargs : dict
        Keywords for Astropy reader.

    """
    from .registries import loader_registry  # Prevent circular import

    logging.info("Attempting to open '{}' using filter '{}'.".format(
            filename, filter))

    name = os.path.basename(filename.name.rstrip(os.sep)).rsplit('.', 1)[0]
    hdulist = fits.open(filename, **kwargs)
    ref = loader_registry.get(filter)

    meta = ref.meta
    header = dict(hdulist[ref.wcs['hdu']].header)
    meta['header'] = header
    wcs = WCS(hdulist[ref.wcs['hdu']].header)

    # Usually, all the data should be in this table
    tab = _read_table(hdulist[ref.data['hdu']], col_idx=ref.data['col'])

    # Read flux column
    data, unit, mask = _read_table_column(tab, ref.data['col'])

    # Find flux unit, if not in column
    if unit is None:
        # Get flux unit from YAML
        if ref.data.get('unit') is not None:
            unit = u.Unit(ref.data['unit'])
        # Get flux unit from header
        else:
            unit = _flux_unit_from_header(meta['header'])

    # Get data mask, if not in column.
    # 0/False = good data (unlike Layers)
    if mask is None:
        mask = np.zeros(data.shape, dtype=np.bool)
    else:
        mask = mask.astype(np.bool)

    # Read in DQ column if it exists
    # 0/False = good (everything else bad)
    if hasattr(ref, 'mask') and ref.mask.get('hdu') is not None:
        if ref.mask['hdu'] == ref.data['hdu']:
            dqtab = tab
        else:
            dqtab = _read_table(
                hdulist[ref.mask['hdu']], col_idx=ref.mask['col'])

        mask2 = _read_table_column(dqtab, ref.mask['col'])[0]  # Data only
        mask |= mask2.astype(np.bool) # Combine with existing mask

    # Wavelength constructed from WCS by default
    dispersion = None
    disp_unit = None

    # Read in wavelength column if it exists
    if hasattr(ref, 'dispersion'):
        if ref.dispersion.get('hdu') is not None:
            if ref.dispersion['hdu'] == ref.data['hdu']:
                wavtab = tab
            else:
                wavtab = _read_table(hdulist[ref.dispersion['hdu']],
                                     col_idx=ref.dispersion['col'])
            dispersion, disp_unit = _read_table_column(
                wavtab, ref.dispersion['col'])[:2]  # Ignore mask

        # Overrides wavelength unit from YAML
        if ref.dispersion.get('unit') is not None:
            disp_unit = u.Unit(ref.dispersion['unit'])

        # If no unit, try to use WCS
        if disp_unit == u.dimensionless_unscaled:
            disp_unit = None

    # Read flux uncertainty
    if hasattr(ref, 'uncertainty') and ref.uncertainty.get('hdu') is not None:
        if ref.uncertainty['hdu'] == ref.data['hdu']:
            errtab = tab
        else:
            errtab = _read_table(
                hdulist[ref.uncertainty['hdu']], col_idx=ref.uncertainty['col'])

        uncertainty = _read_table_column(
            errtab, ref.uncertainty['col'], to_unit=unit)[0]  # Data only
        uncertainty_type = ref.uncertainty.get('type', 'std')
    else:
        uncertainty = np.zeros(data.shape)
        uncertainty_type = 'std'

    # This is dictated by the type of the uncertainty.
    uncertainty = _set_uncertainty(uncertainty, uncertainty_type)

    hdulist.close()

    return Data(name=name, data=data, unit=unit, uncertainty=uncertainty,
                mask=mask, wcs=wcs, dispersion=dispersion,
                dispersion_unit=disp_unit)


# NOTE: This is used by both FITS and ASCII.
def _set_uncertainty(err_array, err_type):
    """Uncertainty is dictated by its type.

    Parameters
    ----------
    err_array : array
        Uncertainty values.

    err_type : {'ivar', 'std'}
        Inverse variance or standard deviation.

    Returns
    -------
    uncertainty : `~astropy.nddata.nduncertainty.StdDevUncertainty`
        Standard deviation uncertainty.

    """
    if err_type == 'ivar':
        err = np.sqrt(1.0 / err_array)
    else:  # 'std'
        err = err_array

    return StdDevUncertainty(err)


def _flux_unit_from_header(header, key='BUNIT'):
    """Get flux unit from header.

    Parameters
    ----------
    header : dict
        Extracted header.

    key : str
        Keyword name.

    Returns
    -------
    unit : `~astropy.units.core.Unit`
        Flux unit. This falls back to default flux unit if look-up failed.

    """
    unitname = header.get(key, default_fluxunit.to_string()).lower()

    # TODO: A more elegant way is to use astropy.units.def_unit()
    if unitname == 'electrons/s':
        unitname = 'electron/s'

    try:
        unit = u.Unit(unitname)
    except ValueError as e:
        unit = default_fluxunit
        logging.warning(str(e))

    return unit


def _read_table(hdu, col_idx=0):
    """Parse FITS table using Astropy first, but use brute force
    if the former fails.

    Astropy parsing is very good at extracting unit and mask, along
    with the data, if FITS table is formatted properly.
    Brute force guarantees the data but provides no unit nor mask.

    Parameters
    ----------
    hdu : obj
        HDU object.

    col_idx : int
        Column index to extract the data directly from HDU.
        This is only used if Astropy parsing fails.

    Returns
    -------
    tab : `~astropy.table.table.Table`
        Parsed table.

    """
    # Let Astropy parse the table for us.
    try:
        tab = Table.read(hdu, format='fits')
    # Build manually if we have to.
    except:
        tab = Table([hdu.data[col_idx].flatten()])

    return tab


def _read_table_column(tab, col_idx, to_unit=None, equivalencies=[]):
    """Read a given Astropy Table column.

    Parameters
    ----------
    tab : `~astropy.table.Table`
        FITS table parsed by Astropy.

    col_idx : int
        Column index.

    to_unit : `~astropy.units.core.Unit` or `None`
        If given, convert data to this unit.

    equivalencies : list
        Astropy unit conversion equivalencies, if needed.
        This might be needed for some flux or wavelength conversions.

    Returns
    -------
    data : array
        1D array of the values.

    unit : `~astropy.units.core.Unit`
        Unit, if any.

    mask : array or `None`
        Mask of the data, if any.

    """
    cols = tab.colnames

    # Special handling for single-column table generated by brute force
    if len(cols) == 1:
        col_idx = 0

    coldat = tab[cols[col_idx]]
    data = coldat.data
    unit = coldat.unit

    # Sometimes, Astropy returns masked column.
    if hasattr(data, 'mask'):
        mask = data.mask.flatten()
        data = data.data.flatten()
    else:
        mask = None
        data = data.flatten()

    # If data has no unit, just assume it is the output unit.
    # Otherwise, perform unit conversion.
    if isinstance(to_unit, u.Unit) and to_unit != u.dimensionless_unscaled:
        unit = to_unit
        if unit != u.dimensionless_unscaled and unit != to_unit:
            data = coldat.to(to_unit, equivalencies).value

    return data, unit, mask


def fits_identify(origin, *args, **kwargs):
    """Check whether given filename is FITS.
    This is used for Astropy I/O Registry.

    """
    return (isinstance(args[0], str) and
            args[0].lower().split('.')[-1] in ['fits', 'fit'])


def ascii_reader(filename, filter, **kwargs):
    """Like :func:`fits_reader` but for ASCII file."""
    from .registries import loader_registry  # Prevent circular import

    name = os.path.basename(filename.name.rstrip(os.sep)).rsplit('.', 1)[0]
    tab = ascii.read(filename, **kwargs)
    cols = tab.colnames
    ref = loader_registry.get(filter)

    meta = ref.meta
    meta['header'] = {}

    # Only loads KEY=VAL comment entries into header
    if 'comments' in tab.meta:
        for s in tab.meta['comments']:
            if '=' not in s:
                continue
            s2 = s.split('=')
            meta['header'][s2[0]] = s2[1]

    wcs = None
    wave = tab[cols[ref.dispersion['col']]]
    dispersion = wave.data
    flux = tab[cols[ref.data['col']]]
    data = flux.data
    uncertainty = np.zeros(data.shape)
    uncertainty_type = 'std'

    if flux.unit is None:
        unit = u.Unit(ref.data.get('unit', default_fluxunit))
    else:
        unit = flux.unit

    if wave.unit is None:
        disp_unit = u.Unit(ref.dispersion.get('unit', default_waveunit))
    else:
        disp_unit = wave.unit

    # 0/False = good data (unlike Layers)
    mask = np.zeros(data.shape, dtype=np.bool)

    if hasattr(ref, 'uncertainty') and ref.uncertainty.get('col') is not None:
        try:
            uncertainty = tab[cols[ref.uncertainty['col']]].data
        except IndexError:
            pass  # Input has no uncertainty column
        else:
            uncertainty_type = ref.uncertainty.get('type', 'std')

    # This is dictated by the type of the uncertainty.
    uncertainty = _set_uncertainty(uncertainty, uncertainty_type)

    if hasattr(ref, 'mask') and ref.mask.get('col') is not None:
        try:
            mask = tab[cols[ref.mask['col']]].data.astype(np.bool)
        except IndexError:
            pass  # Input has no mask column

    return Data(name=str(name), data=data, dispersion=dispersion,
                uncertainty=uncertainty, mask=mask, wcs=wcs,
                unit=unit, dispersion_unit=disp_unit)


def ascii_identify(origin, *args, **kwargs):
    """Check whether given filename is ASCII.
    This is used for Astropy I/O Registry.

    """
    return (isinstance(args[0], str) and
            args[0].lower().split('.')[-1] in ['txt', 'dat'])


def linelist_reader(filename, filter, **kwargs):
    from .registries import loader_registry  # Prevent circular import

    ref = loader_registry.get(filter)

    names_list = []
    start_list = []
    end_list = []
    for k in range(len((ref.columns))):
        name = ref.columns[k]['name']
        names_list.append(name)
        start = ref.columns[k]['start']
        end = ref.columns[k]['end']
        start_list.append(start)
        end_list.append(end)

    tab = ascii.read(filename, format = ref.format,
                     names = names_list,
                     col_starts = start_list,
                     col_ends = end_list)

    return LineList(tab)


def linelist_identify(origin, *args, **kwargs):
    """Check whether given filename is a line list.
    """
    return (isinstance(args[0], str) and
            args[0].lower().split('.')[-1] in ['txt', 'dat'])


# NOTE: Need it this way to prevent circular import.
def register_loaders():
    """Add IO reader/identifier to io registry."""
    from .registries import io_registry

    # FITS
    io_registry.register_reader('fits', Data, fits_reader)
    io_registry.register_identifier('fits', Data, fits_identify)

    # ASCII
    io_registry.register_reader('ascii', Data, ascii_reader)
    io_registry.register_identifier('ascii', Data, ascii_identify)

    # line list
    io_registry.register_reader('ascii', LineList, linelist_reader)
    io_registry.register_identifier('ascii', LineList, linelist_identify)
