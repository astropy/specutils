"""
The core specutils data objects package. This contains the
`~astropy.nddata.NDData`-inherited classes used for storing the spectrum data.
"""
from __future__ import absolute_import

from astropy.nddata import NDIOMixin

from .spectrum1d import *  # noqa
from .spectral_region import *  # noqa
from .spectrum_collection import *  #noqa


class SpectrumList(list, NDIOMixin):
    """
    A list that is used to hold a list of Spectrum1D objects

    The primary purpose of this class is to allow loaders to return a list of
    heterogenous spectra that do have a spectral axis of the same length.
    """
