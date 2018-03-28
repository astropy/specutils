from __future__ import absolute_import, division

import logging
from functools import wraps

from astropy.io import registry as io_registry

from ..spectra.spectrum1d import Spectrum1D


def data_loader(label, identifier=None, dtype=Spectrum1D):
    """
    Wraps a function that can be added to an `~astropy.io.registry` for custom
    file reading.

    Parameters
    ----------
    label : str
        The label given to the function inside the registry.
    identifier : func
        The identified function used to verify that a file is to use a
        particular file.
    dtype : class
        A class reference for which the data loader should be store.
    """
    def decorator(func):
        io_registry.register_reader(label, dtype, func)
        io_registry.register_identifier(label, dtype, identifier)

        logging.debug("Successfully loaded reader \"{}\".".format(label))

        @wraps(func)
        def wrapper(*args, **kwargs):
            return func(*args, **kwargs)
        return wrapper
    return decorator


def custom_writer(label, dtype=Spectrum1D):
    def decorator(func):
        io_registry.register_writer(label, Spectrum1D, func)

        @wraps(func)
        def wrapper(*args, **kwargs):
            return func(*args, **kwargs)
        return wrapper
    return decorator
