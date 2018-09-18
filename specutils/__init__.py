# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Specutils: an astropy package for spectroscopy.
"""

# Packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *
# ----------------------------------------------------------------------------

# Enforce Python version check during package import.
# This is the same check as the one at the top of setup.py
import sys

__minimum_python_version__ = "3.5"

class UnsupportedPythonError(Exception):
    pass

if sys.version_info < tuple((int(val) for val in __minimum_python_version__.split('.'))):
    raise UnsupportedPythonError("packagename does not support Python < {}".format(__minimum_python_version__))

if not _ASTROPY_SETUP_:
    # For egg_info test builds to pass, put package imports here.

    # Allow loading spectrum object from top level module
    from .spectra import *

    # Load the default IO functions
    from .io.default_loaders import *  # noqa


    def load_user():
        import os
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

    load_user()
