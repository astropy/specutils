from ..wcs_adapter import WCSAdapter, WCSAxes
from astropy.modeling import models
from gwcs.wcs import WCS

__all__ = ['GWCSAdapter']


class GWCSAdapter(WCSAdapter):
    wrapped_class = WCS
    axes = None

    def __init__(self, wcs):
        super(GWCSAdapter, self).__init__(wcs)

        self.axes = WCSAxes(
            spectral=self.wcs.forward_transform | models.Mapping((2,))
        )

    def world_to_pix(self):
        pass

    def pix_to_world(self):
        pass