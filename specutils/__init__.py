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

    # Load the IO functions
    from .io.default_loaders import *  # noqa
    from .io.registers import _load_user_io
    _load_user_io()

__citation__ = 'https://doi.org/10.5281/zenodo.1421356'
