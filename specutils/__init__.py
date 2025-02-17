# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Specutils: an astropy package for spectroscopy.
"""
import os

from astropy import config as _config

__citation__ = 'https://doi.org/10.5281/zenodo.1421356'

try:
    from .version import version as __version__
except ImportError:
    __version__ = ''

# Allow loading spectrum object from top level module
from .spectra import *  # noqa

# Load the IO functions
from .io.default_loaders import *  # noqa
from .io.registers import _load_user_io  # noqa
_load_user_io()


class Conf(_config.ConfigNamespace):
    """Configuration parameters for specutils."""

    do_continuum_function_check = _config.ConfigItem(
        True,
        'Whether to check the spectrum baseline value is close'
        'to zero. If it is not within ``threshold`` then a warning is raised.'
    )


conf = Conf()

# Clean up namespace
del os
del _config
del _load_user_io
