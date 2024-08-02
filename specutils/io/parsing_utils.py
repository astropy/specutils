import numpy as np
import os
import re
import urllib
import io
import contextlib

from astropy.io import fits
from astropy.nddata import StdDevUncertainty
from astropy.utils.exceptions import AstropyUserWarning
import astropy.units as u
import warnings

from specutils.spectra import Spectrum1D


@contextlib.contextmanager
def read_fileobj_or_hdulist(*args, **kwargs):
    """ Context manager for reading a filename or file object

    Returns
    -------
    hdulist : :class:`~astropy.io.fits.HDUList`
        Provides a generator-iterator representing the open file object handle.
    """
    # Access the fileobj or filename arg
    # Do this so identify functions are useable outside of Spectrum1d.read context
    try:
        fileobj = args[2]
    except IndexError:
        fileobj = args[0]

    if isinstance(fileobj, fits.hdu.hdulist.HDUList):
        if fits.util.fileobj_closed(fileobj):
            hdulist = fits.open(fileobj.name, **kwargs)
        else:
            hdulist = fileobj
    elif isinstance(fileobj, io.BufferedReader):
        hdulist = fits.open(fileobj)
    else:
        hdulist = fits.open(fileobj, **kwargs)

    try:
        yield hdulist

    # Cleanup even after identifier function has thrown an exception: rewind generic file handles.
    finally:
        if not isinstance(fileobj, fits.hdu.hdulist.HDUList):
            try:
                fileobj.seek(0)
            except (AttributeError, io.UnsupportedOperation):
                hdulist.close()


def spectrum_from_column_mapping(table, column_mapping, wcs=None, verbose=False):
    """
    Given a table and a mapping of the table column names to attributes
    on the Spectrum1D object, parse the information into a Spectrum1D.

    Parameters
    ----------
    table : :class:`~astropy.table.Table`
        The table object (e.g. returned from ``Table.read('data_file')``).

    column_mapping : dict
        A dictionary describing the relation between the table columns
        and the arguments of the `Spectrum1D` class, along with unit
        information. The dictionary keys should be the table column names
        while the values should be a two-tuple where the first element is the
        associated `Spectrum1D` keyword argument, and the second element is the
        unit for the file column (or ``None`` to take unit from the table header)::

            column_mapping = {'FLUX': ('flux', 'Jy'),
                              'WAVE': ('spectral_axis', 'um')}

    wcs : :class:`~astropy.wcs.WCS` or :class:`gwcs.WCS`
        WCS object passed to the Spectrum1D initializer.

    verbose : bool
        Print extra info.

    Returns
    -------
    :class:`~specutils.Spectrum1D`
        The spectrum with 'spectral_axis', 'flux' and optionally 'uncertainty'
        as identified by `column_mapping`.
    """
    spec_kwargs = {}

    # Associate columns of the file with the appropriate spectrum1d arguments
    for col_name, (kwarg_name, cm_unit) in column_mapping.items():
        # If the table object couldn't parse any unit information,
        # fallback to the column mapper defined unit
        tab_unit = table[col_name].unit

        if tab_unit and cm_unit is not None:
            # If the table unit is defined, retrieve the quantity array for
            # the column
            kwarg_val = u.Quantity(table[col_name], tab_unit)

            # Attempt to convert the table unit to the user-defined unit.
            if verbose:
                print(f"Attempting auto-convert of table unit '{tab_unit}' to "
                      f"user-provided unit '{cm_unit}'.")

            if not isinstance(cm_unit, u.Unit):
                cm_unit = u.Unit(cm_unit)
            if cm_unit.physical_type in ('length', 'frequency', 'energy'):
                # Spectral axis column information
                kwarg_val = kwarg_val.to(cm_unit, equivalencies=u.spectral())
            elif 'spectral flux' in str(cm_unit.physical_type):
                # Flux/error column information
                kwarg_val = kwarg_val.to(cm_unit, equivalencies=u.spectral_density(1 * u.AA))
        elif tab_unit:
            # The user has provided no unit in the column mapping, so we
            # use the unit as defined in the table object.
            kwarg_val = u.Quantity(table[col_name], tab_unit)
        elif cm_unit is not None:
            # In this case, the user has defined a unit in the column mapping
            # but no unit has been defined in the table object.
            kwarg_val = u.Quantity(table[col_name], cm_unit)
        else:
            # Neither the column mapping nor the table contain unit information.
            # This may be desired e.g. for the mask or bit flag arrays.
            kwarg_val = table[col_name]

        # Transpose > 1D data to row-major format
        if kwarg_val.ndim > 1:
            kwarg_val = kwarg_val.T

        spec_kwargs.setdefault(kwarg_name, kwarg_val)

    # Ensure that the uncertainties are a subclass of NDUncertainty
    if spec_kwargs.get('uncertainty') is not None:
        spec_kwargs['uncertainty'] = StdDevUncertainty(
            spec_kwargs.get('uncertainty'))

    return Spectrum1D(**spec_kwargs, wcs=wcs, meta={'header': table.meta})


def generic_spectrum_from_table(table, wcs=None):
    """
    Load spectrum from an Astropy table into a Spectrum1D object.
    Uses the following logic to figure out which column is which:

     * Spectral axis (dispersion) is the first column with units
     compatible with ``u.spectral()`` or with length units such as 'pix'.
     Need not be present, if a valid ``wcs`` parameter is passed.

     * Flux is taken from the first column with units compatible with
     ``u.spectral_density()``, or with other likely culprits such as
     'adu' or 'cts/s'.

     * Uncertainty comes from the next column with the same units as flux.

    Parameters
    ----------
    table : :class:`~astropy.table.Table`
        Table containing a column of ``flux``, and optionally ``spectral_axis``
        and ``uncertainty`` as defined above.
    wcs : :class:`~astropy.wcs.WCS`
        A FITS WCS object. If this is present, the machinery will fall back
        and default to using the ``wcs`` to find the dispersion information.

    Returns
    -------
    :class:`~specutils.Spectrum1D`
        The spectrum that is represented by the data from the columns
        as automatically identified above.

    Raises
    ------
    Warns if uncertainty has zeros or negative numbers.
    Raises IOError if it can't figure out the columns.

    """
    # Local function to find the wavelength or frequency column
    def _find_spectral_axis_column(table, columns_to_search):
        """
        Figure out which column in a table holds the spectral axis (dispersion).
        Take the first column that has units compatible with u.spectral()
        equivalencies. If none meet that criterion, look for other likely
        length units such as 'pix'.
        """
        additional_valid_units = [u.Unit('pix')]
        found_column = None

        # First, search for a column with units compatible with Angstroms
        for c in columns_to_search:
            try:
                table[c].to("AA", equivalencies=u.spectral())
                found_column = c
                break
            except Exception:
                continue

        # If no success there, check for other possible length units
        if found_column is None:
            for c in columns_to_search:
                if table[c].unit in additional_valid_units:
                    found_column = c
                    break

        return found_column

    # Local function to find the flux column
    def _find_spectral_column(table, columns_to_search, spectral_axis):
        """
        Figure out which column in a table holds the fluxes or uncertainties.
        Take the first column that has units compatible with
        u.spectral_density() equivalencies. If none meet that criterion,
        look for other likely length units such as 'adu' or 'cts/s'.
        """
        additional_valid_units = [u.Unit('adu'), u.Unit('ct/s'), u.Unit('count')]
        found_column = None

        # First, search for a column with units compatible with Jansky
        for c in columns_to_search:
            try:
                # Check for multi-D flux columns
                if table[c].ndim == 1:
                    spec_ax = spectral_axis
                else:
                    # Assume leading dimension corresponds to spectral_axis
                    spec_shape = np.ones(table[c].ndim, dtype=int)
                    spec_shape[0] = -1
                    spec_ax = spectral_axis.reshape(spec_shape)
                table[c].to("Jy", equivalencies=u.spectral_density(spec_ax))
                found_column = c
                break
            except Exception:
                continue

        # If no success there, check for other possible flux units
        if found_column is None:
            for c in columns_to_search:
                if table[c].unit in additional_valid_units:
                    found_column = c
                    break

        return found_column

    # Make a copy of the column names so we can remove them as they are found
    colnames = table.colnames.copy()

    # Use the first column that has spectral unit as the dispersion axis
    spectral_axis_column = _find_spectral_axis_column(table, colnames)

    if spectral_axis_column is None and wcs is None:
        raise IOError("Could not identify column containing the wavelength, frequency or energy")
    elif wcs is not None:
        spectral_axis = None
    else:
        spectral_axis = table[spectral_axis_column].to(table[spectral_axis_column].unit)
        colnames.remove(spectral_axis_column)

    # Use the first column that has a spectral_density equivalence as the flux
    flux_column = _find_spectral_column(table, colnames, spectral_axis)
    if flux_column is None:
        raise IOError("Could not identify column containing the flux")
    flux = table[flux_column].to(table[flux_column].unit)
    colnames.remove(flux_column)
    # For > 1D data transpose to row-major format
    if flux.ndim > 1:
        flux = flux.T

    # Use the next column with the same units as flux as the uncertainty
    # Interpret it as a standard deviation and check if it has zeros or negative values
    err_column = None
    for c in colnames:
        if table[c].unit == table[flux_column].unit:
            err_column = c
            break
    if err_column is not None:
        if table[err_column].ndim > 1:
            err = table[err_column].T
        elif flux.ndim > 1:  # Repeat uncertainties over all flux columns
            err = np.tile(table[err_column], flux.shape[0], 1)
        else:
            err = table[err_column]
        err = StdDevUncertainty(err.to(err.unit))
        if np.min(table[err_column]) <= 0.:
            warnings.warn("Standard Deviation has values of 0 or less", AstropyUserWarning)
    else:
        err = None

    # Check for mask
    if 'mask' in table.colnames:
        mask = table['mask']
        if mask.ndim > 1:
            mask = mask.T
    else:
        mask = None

    # Create the Spectrum1D object and return it
    if wcs is not None or spectral_axis_column is not None and flux_column is not None:
        # For > 1D spectral axis transpose to row-major format and return SpectrumCollection
        spectrum = Spectrum1D(flux=flux, spectral_axis=spectral_axis,
                              uncertainty=err, meta={'header': table.meta}, wcs=wcs,
                              mask=mask)

    return spectrum


def _fits_identify_by_name(origin, fileinp, *args,
                           pattern=r'(?i).*\.fit[s]?$', **kwargs):
    """
    Check whether input file is FITS and matches a given name pattern.
    Utility function to construct an `identifier` for Astropy I/O Registry.

    Parameters
    ----------
    fileinp : str or file-like object
        FITS file name or object (provided from name by Astropy I/O Registry).
    pattern : regex str or re.Pattern
        File name pattern to be matched.
        Note: loaders should define a pattern sufficiently specific for their
        spectrum file types to avoid ambiguous/multiple matches.
    """
    fileobj = None
    filepath = None
    if pattern is None:
        pattern = r''
    _spec_pattern = re.compile(pattern)

    if isinstance(fileinp, str):
        filepath = fileinp
        try:
            fileobj = open(filepath, mode='rb')
        except FileNotFoundError:
            # Check if path points to valid url
            try:
                fileinp = urllib.request.urlopen(filepath)
            except ValueError:
                return False
    elif fits.util.isfile(fileinp):
        fileobj = fileinp
        filepath = fileobj.name

    # Check for `urlopen` object - can only probe content if seekable
    if hasattr(fileinp, 'url') and hasattr(fileinp, 'seekable'):
        filepath = urllib.parse.unquote(fileinp.url)
        if fileinp.seekable():
            fileobj = fileinp

    check = (_spec_pattern.match(os.path.basename(filepath)) is not None and
             fits.connect.is_fits(origin, filepath, fileobj, *args))

    if fileobj is not None:
        fileobj.close()

    return check
