# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Specutils: an astropy package for spectroscopy.
"""

# Packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *   # noqa
from astropy import config as _config
# ----------------------------------------------------------------------------

if not _ASTROPY_SETUP_:  # noqa
    # For egg_info test builds to pass, put package imports here.

    # Allow loading spectrum object from top level module
    from .spectra import *  # noqa

    # Load the IO functions
    from .io.default_loaders import *  # noqa
    from .io.registers import _load_user_io
    _load_user_io()

__citation__ = 'https://doi.org/10.5281/zenodo.1421356'


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for specutils.
    """

    do_continuum_function_check = _config.ConfigItem(
        True,
        'Whether to check the spectrum baseline value is close'
        'to zero. If it is not within ``threshold`` then a warning is raised.'
    )


conf = Conf()
