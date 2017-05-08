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
            longitude=self.wcs.sub([WCSSUB_LONGITUDE]),
            latitude=self.wcs.sub([WCSSUB_LATITUDE]),
            cubeface=self.wcs.sub([WCSSUB_CUBEFACE]),
            spectral=self.wcs.sub([WCSSUB_SPECTRAL]),
            stokes=self.wcs.sub([WCSSUB_STOKES]),
            celestial=self.wcs.sub([WCSSUB_CELESTIAL])
        )

    def world_to_pix(self):
        pass

    def pix_to_world(self):
        pass