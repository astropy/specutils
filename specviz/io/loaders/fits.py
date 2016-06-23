from specviz.interfaces.decorators import data_loader

from astropy import units as u
import numpy as np
from astropy import units as u
from astropy.io import ascii, fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.nddata import StdDevUncertainty

from specutils import Spectrum1D

import logging
import os

# Loader automatically falls back to these units for some cases
default_waveunit = u.Unit('Angstrom')
default_fluxunit = u.Unit('erg / (Angstrom cm2 s)')


def fits_identify(origin, *args, **kwargs):
    """Check whether given filename is FITS.
    This is used for Astropy I/O Registry.

    """
    return (isinstance(args[0], str) and
            args[0].lower().split('.')[-1] in ['fits', 'fit'])


@data_loader(label="fits", identifier=fits_identify)
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
    logging.info("Attempting to open '{}' using filter '{}'.".format(
        filename, filter))

    name = os.path.basename(filename.name.rstrip(os.sep)).rsplit('.', 1)[0]
    hdulist = fits.open(filename, **kwargs)

    wcs = WCS(hdulist[0].header)

    # Use Astropy Tables to search the fits file for real data
    t = Table.read(name)

    if len(t.colnames) > 1:
        data_col = [x for x in t.colnames if x.lower() in ['flux', 'data']]
        data_col = data_col[0] if len(data_col) > 0 else 1
        data = t[data_col]

        wave_col = [x for x in t.colenames if x.lower() in ['wave',
                                                            'wavelength']]


    hdulist.close()

    return Spectrum1D.from_array(name=name, data=data, unit=unit,
                                 uncertainty=uncertainty, mask=mask,
                                 wcs=wcs, dispersion=dispersion,
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