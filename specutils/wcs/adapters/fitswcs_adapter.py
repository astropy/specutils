from ..wcs_adapter import WCSAdapter, WCSAxes

from astropy.wcs import WCS
from astropy.wcs import WCSSUB_SPECTRAL, WCSSUB_LONGITUDE, WCSSUB_LATITUDE, \
    WCSSUB_CUBEFACE, WCSSUB_STOKES, WCSSUB_CELESTIAL


__all__ = ['FITSWCSAdapter']


class FITSWCSAdapter(WCSAdapter):
    wrapped_class = WCS
    axes = None

    def __init__(self, wcs):
        super(FITSWCSAdapter, self).__init__(wcs)

        self.axes = WCSAxes(
            longitude=FITSWCSAdapter(self.wcs.sub([WCSSUB_LONGITUDE])),
            latitude=FITSWCSAdapter(self.wcs.sub([WCSSUB_LATITUDE])),
            cubeface=FITSWCSAdapter(self.wcs.sub([WCSSUB_CUBEFACE])),
            spectral=FITSWCSAdapter(self.wcs.sub([WCSSUB_SPECTRAL])),
            stokes=FITSWCSAdapter(self.wcs.sub([WCSSUB_STOKES])),
            celestial=FITSWCSAdapter(self.wcs.sub([WCSSUB_CELESTIAL]))
        )

    def world_to_pix(self):
        pass

    def pix_to_world(self):
        pass

    def bin_edges(self):
        # the WCS doesn't know about its own pixel array
        edge_indices = list(self.wcs.pixel_indices-0.5) + \
                       [self.wcs.pixel_indices[-1]+0.5]

        return self.pix2world(edge_indices, 0)