# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
The specutils package implements base classes and utilities for
interacting with astronomical spectra in Python and the Astropy
project.  It is intended for eventual merger with the `astropy`
package, but for now is being developed independently.
"""

# Affiliated packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *
# ----------------------------------------------------------------------------

if not _ASTROPY_SETUP_:

    from .spectrum1d import Spectrum1D
    from . import io
