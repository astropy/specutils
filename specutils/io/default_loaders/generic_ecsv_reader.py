import logging
import six
import os
from ...spectra import Spectrum1D
from ..registers import data_loader
from astropy.table import Table
from ..parsing_utils import generic_spectrum_from_table


def identify_ecsv(origin, *args, **kwargs):
    """Check if it's an ECSV file."""
    return (isinstance(args[0], six.string_types) and
            os.path.splitext(args[0].lower())[1] == '.ecsv')


@data_loader("generic-ecsv", identifier=identify_ecsv, dtype=Spectrum1D)
def generic_ecsv(file_name, **kwargs):
    """
    Read a spectrum from an ECSV file, using generic_spectrum_from_table_loader()
    to try to figure out which column is which.

    Parameters
    ----------
    file_name: str
        The path to the ECSV file

    Returns
    -------
    data: Spectrum1D
        The spectrum that is represented by the data in this table.
    """
    logging.info("Spectrum file looks like 'generic-ecsv'.")

    table = Table.read(file_name, format='ascii.ecsv')

    spectrum = generic_spectrum_from_table(table, **kwargs)

    return spectrum
