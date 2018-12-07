import numpy as np
from astropy.table import Table
import astropy.units as u
from astropy.nddata import StdDevUncertainty
from astropy.utils.exceptions import AstropyUserWarning
import warnings
import logging

from specutils.spectra import Spectrum1D


def spectrum_from_column_mapping(table, column_mapping, wcs=None):
    """
    Given a table and a mapping of the table column names to attributes
    on the Spectrum1D object, parse the information into a Spectrum1D.

    Parameters
    ----------
    table : :class:`~astropy.table.Table`
        The table object returned from parsing the data file.
    column_mapping : dict
        A dictionary describing the relation between the file columns
        and the arguments of the `Spectrum1D` class, along with unit
        information. The dictionary keys should be the file column names
        while the values should be a two-tuple where the first element is the
        associated `Spectrum1D` keyword argument, and the second element is the
        unit for the file column::

            column_mapping = {'FLUX': ('flux', 'Jy')}

    wcs : :class:`~astropy.wcs.WCS` or :class:`gwcs.WCS`
        WCS object passed to the Spectrum1D initializer.
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
            logging.debug("Attempting auto-convert of table unit '%s' to "
                          "user-provided unit '%s'.", tab_unit, cm_unit)

            if cm_unit.physical_type in ('length', 'frequency'):
                # Spectral axis column information
                kwarg_val = kwarg_val.to(cm_unit, equivalence=u.spectral())
            elif 'spectral flux' in cm_unit.physical_type:
                # Flux/error column information
                kwarg_val = kwarg_val.to(
                    cm_unit, equivalencies=u.spectral_density(1 * u.AA))
        elif cm_unit is not None:
            # In this case, the user has defined a unit in the column mapping
            # but no unit has been defined in the table object.
            kwarg_val = u.Quantity(table[col_name], cm_unit)
        else:
            # Neither the column mapping nor the table contain unit information.
            # This may be desired e.g. for the mask or bit flag arrays.
            kwarg_val = table[col_name]

        spec_kwargs.setdefault(kwarg_name, kwarg_val)

    # Ensure that the uncertainties are a subclass of NDUncertainty
    if spec_kwargs.get('uncertainty') is not None:
        spec_kwargs['uncertainty'] = StdDevUncertainty(
            spec_kwargs.get('uncertainty'))

    return Spectrum1D(**spec_kwargs, wcs=wcs, meta=table.meta)


def generic_spectrum_from_table(table, wcs=None, **kwargs):
    """
    Load spectrum from an Astropy table into a Spectrum1D object.
    Uses the following logic to figure out which column is which:

     * Spectral axis (dispersion) is the first column with units
     compatible with u.spectral() or with length units such as 'pix'.

     * Flux is taken from the first column with units compatible with
     u.spectral_density(), or with other likely culprits such as
     'adu' or 'cts/s'.

     * Uncertainty comes from the next column with the same units as flux.

    Parameters
    ----------
    file_name: str
        The path to the ECSV file
    wcs : :class:`~astropy.wcs.WCS`
        A FITS WCS object. If this is present, the machinery will fall back
        to using the wcs to find the dispersion information.

    Returns
    -------
    data: Spectrum1D
        The spectrum that is represented by the data in this table.

    Raises
    ------
    Warns if uncertainty has zeros or negative numbers.
    Raises IOError if it can't figure out the columns.

    """
    # Local function to find the wavelength or frequency column
    def _find_spectral_axis_column(table,columns_to_search):
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
                table[c].to("AA",equivalencies=u.spectral())
                found_column = c
                break
            except:
                continue

        # If no success there, check for other possible length units
        if found_column is None:
            for c in columns_to_search:
                if table[c].unit in additional_valid_units:
                    found_column = c
                    break

        return found_column

    # Local function to find the flux column
    def _find_spectral_column(table,columns_to_search,spectral_axis):
        """
        Figure out which column in a table holds the fluxes or uncertainties.
        Take the first column that has units compatible with
        u.spectral_density() equivalencies. If none meet that criterion,
        look for other likely length units such as 'adu' or 'cts/s'.
        """
        additional_valid_units = [u.Unit('adu'),u.Unit('ct/s')]
        found_column = None

        # First, search for a column with units compatible with Janskies
        for c in columns_to_search:
            try:
                table[c].to("Jy",equivalencies=u.spectral_density(spectral_axis))
                found_column = c
                break
            except:
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
    flux_column = _find_spectral_column(table,colnames,spectral_axis)
    if flux_column is None:
        raise IOError("Could not identify column containing the flux")
    flux = table[flux_column].to(table[flux_column].unit)
    colnames.remove(flux_column)

    # Use the next column with the same units as flux as the uncertainty
    # Interpret it as a standard deviation and check if it has zeros or negative values
    err_column = None
    for c in colnames:
        if table[c].unit == table[flux_column].unit:
            err_column = c
            break
    if err_column is not None:
        err = StdDevUncertainty(table[err_column].to(table[err_column].unit))
        if np.min(table[err_column]) <= 0.:
            warnings.warn("Standard Deviation has values of 0 or less", AstropyUserWarning)

    # Create the Spectrum1D object and return it
    if wcs is not None or spectral_axis_column is not None and flux_column is not None:
        if err_column is not None:
            spectrum = Spectrum1D(flux=flux, spectral_axis=spectral_axis,
                                  uncertainty=err, meta=table.meta, wcs=wcs)
        else:
            spectrum = Spectrum1D(flux=flux, spectral_axis=spectral_axis,
                                  meta=table.meta, wcs=wcs)

    return spectrum
