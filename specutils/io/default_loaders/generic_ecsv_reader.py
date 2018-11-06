import logging
import six
import os
from ...spectra import Spectrum1D
from ..registers import data_loader
from astropy.table import Table
from ..generic_spectrum_from_table import generic_spectrum_from_table

def identify_ecsv(origin, *args, **kwargs):
# check if file can be opened with this reader
# args[0] = filename
    return (isinstance(args[0], six.string_types) and
    # check if file is .ecsv (would be nice to validate the header too)
        os.path.splitext(args[0].lower())[1] == '.ecsv'
        )

@data_loader("generic-ecsv", identifier=identify_ecsv, dtype=Spectrum1D)
def generic_ecsv(file_name, **kwargs):
    """
    Read a spectrum from an ECSV file, using generic_spectrum_from_table_loader()
    to try to figure out which column is which.
    """
    logging.info("Spectrum file looks like generic-ecsv")
    table = Table.read(file_name, format='ascii.ecsv')
    spectrum = generic_spectrum_from_table(table, **kwargs)
    return spectrum
