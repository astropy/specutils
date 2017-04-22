import numpy as np
from astropy.modeling import models
from astropy import units as u
from astropy import coordinates as coords
from gwcs import wcs, utils
from gwcs import coordinate_frames as cf
from astropy.modeling.tabular import Tabular1D
import astropy.wcs.WCS

def tabular_wcs(xarray):

    coordinateframe = cf.CoordinateFrame(naxes=1, axes_type=('SPECTRAL',), axes_order=(0,))
    specframe = cf.SpectralFrame(unit=xarray.unit, axes_order=(0,))
    transform = Tabular1D(np.arange(len(xarray)), xarray.value)

    tabular_gwcs = wcs.WCS([(coordinateframe, transform), (specframe, None)])

    return tabular_gwcs

class TabularGWCSWrapper(object):
    def __init__(self, xarray):
        self.xarray = xarray
        self.gwcs = tabular_wcs(xarray)

    def pix_to_world(self, pixel_array):
        return self.gwcs(pixel_array, with_bounding_box=False)

    def world_to_pix(self, world_coordinate_array):
        return self.gwcs.backward_transform(world_coordinate_array)

    def with_new_unit(self, new_unit, equivalencies=None):
        return self.__class__(self.xarray.to(new_unit, equivalencies=equivalencies))

    def bin_edges(self):
        raise NotImplementedError("What interpretation can we assign to tabular bin edges?"
                                  "Maybe we have a different object that understands bin edges"
                                  " and just raise exceptions for tabular ones because we don't "
                                  "understand what the bin edges are by default?")

class FITSWCSWrapper(astropy.wcs.WCS):
    def bin_edges(self):
        # the WCS doesn't know about its own pixel array
        edge_indices = list(self.pixel_indices-0.5) + [self.pixel_indices[-1]+0.5]
        return self.wcs_pix2world(edge_indices, 0)
