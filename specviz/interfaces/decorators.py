from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import logging

import astropy.io.registry as io_registry

from ..core.data import GenericSpectrum1D


def data_loader(label, identifier):
    """
    A decorator that registers a function and identifies with an Astropy io
    registry object.

    Parameters
    ----------
    load_func : function
        Function added to the registry in order to read data files.
    """
    def decorator(reader_func):
        def wrapper():
            logging.info("Added {} to loader registry.".format(label))
            io_registry.register_reader(label, GenericSpectrum1D, reader_func)
            io_registry.register_identifier(label, GenericSpectrum1D, identifier)

            return

        return wrapper

    return decorator