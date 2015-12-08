from __future__ import absolute_import, division, print_function

import os
import numpy as np

from ..core.data import Data
from .registries import loader_factory
from astropy.io import registry, fits
from astropy.wcs import WCS


def fits_reader(filename, filter, **kwargs):
    """
    This generic function will query the loader factory which has already loaded the yaml configuration files in an
    attempt to parse the associated fits file.

    Parameters
    ----------
    filename
    filter
    kwargs

    Returns
    -------

    """
    hdulist = fits.open(filename, **kwargs)
    ref = loader_factory.get(filter)

    wcs = WCS(hdulist[ref.wcs['hdu']].header)
    data = hdulist[ref.data['hdu']].data
    uncertainty = hdulist[ref.uncertainty['hdu']].data
    mask = hdulist[ref.mask['hdu']].data if ref.mask['hdu'] is not None else np.zeros(shape=data.shape)

    return Data(data=)


def fits_identify(origin, *args, **kwargs):
    return isinstance(args[0], basestring) and \
           args[0].lower().split('.')[-1] in ['fits', 'fit']


registry.register_reader('fits', Data, fits_reader)
registry.register_identifier('fits', Data, fits_identify)
