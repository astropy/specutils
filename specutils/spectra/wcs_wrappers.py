import numpy as np
from astropy.modeling import models
from astropy import units as u
from astropy import coordinates as coords
from gwcs import wcs, utils
from gwcs import coordinate_frames as cf
from astropy.modeling.tabular import Tabular1D

def tabular_wcs(xarray):

    coordinateframe = cf.CoordinateFrame(naxes=1, axes_type=('SPECTRAL',), axes_order=(0,))
    specframe = cf.SpectralFrame(unit=xarray.unit, axes_order=(0,))
    transform = Tabular1D(np.arange(len(xarray)), xarray.value)

    tabular_gwcs = wcs.WCS([(coordinateframe, transform), (specframe, None)])

    return tabular_gwcs
