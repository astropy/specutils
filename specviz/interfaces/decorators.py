from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import astropy.io.registry as io_registry


def data_loader(label, identifier, cls):
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
            print("ADDED")
            io_registry.register_reader(label, cls, reader_func)
            io_registry.register_identifier(label, cls, identifier)

            return

        return wrapper

    return decorator