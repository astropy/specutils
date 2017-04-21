# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This is an Astropy affiliated package.
"""

# Affiliated packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *
# ----------------------------------------------------------------------------
import sys
import os
import logging

logging.basicConfig(level=logging.INFO)

# Load the default IO functions
from .io.default_loaders import *

def load_user():
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
            except ModuleNotFoundError:
                pass

load_user()
