"""
A module containing the mechanics of the specutils io registry.
"""
import os
import logging
from functools import wraps

from astropy.io import registry as io_registry

from ..spectra.spectrum1d import Spectrum1D


__all__ = ['data_loader', 'custom_writer']


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

def _load_user_io():
    # Get the path relative to the user's home directory
    path = os.path.expanduser("~/.specutils")

    # If the directory doesn't exist, create it
    if not os.path.exists(path):
        os.mkdir(path)

    # Import all python files from the directory
    for file in os.listdir(path):
        if not file.endswith("py"):
            continue

        try:
            import importlib.util as util

            spec = util.spec_from_file_location(file[:-3],
                                                os.path.join(path, file))
            mod = util.module_from_spec(spec)
            spec.loader.exec_module(mod)
        except ImportError:
            from importlib import import_module

            sys.path.insert(0, path)

            try:
                import_module(file[:-3])
            except ModuleNotFoundError:  # noqa
                pass
