"""
A module containing the mechanics of the specutils io registry.
"""
import os
import logging
from functools import wraps

from astropy.io import registry as io_registry

from ..spectra import Spectrum1D, SpectrumList


__all__ = ['data_loader', 'custom_writer', 'get_loaders_by_extension']


def data_loader(label, identifier=None, dtype=Spectrum1D, extensions=None,
                priority=0):
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
    extensions : list
        A list of file extensions this loader supports loading from. In the
        case that no identifier function is defined, but a list of file
        extensions is, a simple identifier function will be created to check
        for consistency with the extensions.
    priority : int
        Set the priority of the loader. Currently influences the sorting of the
        returned loaders for a dtype.
    """
    def decorator(func):
        io_registry.register_reader(label, dtype, func)

        # If the identifier is not defined, but the extensions are, create
        # a simple identifier based off file extension.
        if identifier is None and extensions is not None:
            io_registry.register_identifier(
                label, dtype, lambda *args, **kwargs: any([args[0].endswith(x)
                                                           for x in extensions]))
        else:
            io_registry.register_identifier(label, dtype, identifier)

        # Include the file extensions as attributes on the function object
        func.extensions = extensions

        # Include priority on the loader function attribute
        func.priority = priority

        # Sort the io_registry based on priority
        sorted_loaders = sorted(io_registry._readers.items(),
                                key=lambda item: getattr(item[1], 'priority', 0))

        # Update the registry with the sorted dictionary
        io_registry._readers.clear()
        io_registry._readers.update(sorted_loaders)

        logging.debug("Successfully loaded reader \"{}\".".format(label))

        # Automatically register a SpectrumList reader for any data_loader that
        # reads Spectrum1D objects. TODO: it's possible that this
        # functionality should be opt-in rather than automatic.
        if dtype is Spectrum1D:
            def load_spectrum_list(*args, **kwargs):
                return SpectrumList([ func(*args, **kwargs) ])

            # Add these attributes to the SpectrumList reader as well
            load_spectrum_list.extensions = extensions
            load_spectrum_list.priority = priority

            io_registry.register_reader(label, SpectrumList, load_spectrum_list)
            logging.debug("Created SpectrumList reader for \"{}\".".format(label))

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


def get_loaders_by_extension(extension):
    """
    Retrieve a list of loader labels associated with a given extension.

    Parameters
    ----------
    extension : str
        The extension for which associated loaders will be matched against.

    Returns
    -------
    loaders : list
        A list of loader names that are associated with the extension.
    """
    return [fmt for (fmt, cls), func in io_registry._readers.items()
            if issubclass(cls, Spectrum1D) and
            func.extensions is not None and
            extension in func.extensions]


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
