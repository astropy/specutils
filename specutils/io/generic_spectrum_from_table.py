import numpy as np
from astropy.table import Table
import astropy.units as u
from astropy.nddata import StdDevUncertainty
from astropy.utils.exceptions import AstropyUserWarning
import warnings

from specutils.spectra import Spectrum1D

def generic_spectrum_from_table(table, **kwargs):
    """
    Load spectrum an Astropy table.
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

    Returns
    -------
    data: Spectrum1D
        The data.

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
        if found_column == None:
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
        if found_column == None:
            for c in columns_to_search:
                if table[c].unit in additional_valid_units:
                    found_column = c
                    break

        return found_column

    # Make a copy of the column names so we can remove them as they are found
    colnames = table.colnames.copy()

    # Use the first column that has spectral unit as the dispersion axis
    spectral_axis_column = _find_spectral_axis_column(table,colnames)
    if spectral_axis_column is None:
       raise IOError("Could not identify column containing the wavelength, frequency or energy")
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
    if spectral_axis_column is not None and flux_column is not None:
       if err_column is not None:
           spectrum = Spectrum1D(flux=flux, spectral_axis=spectral_axis,
               uncertainty=err,meta=table.meta)
       else:
           spectrum = Spectrum1D(flux=flux, spectral_axis=spectral_axis,meta=table.meta)
    return spectrum
