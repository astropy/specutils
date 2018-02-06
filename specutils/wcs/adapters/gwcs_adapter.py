from ..wcs_adapter import WCSAdapter
from gwcs.wcs import WCS

__all__ = ['GWCSAdapter']


class GWCSAdapter(WCSAdapter):
    """
    Adapter class that adds support for GWCS objects.
    """
    wrapped_class = WCS
    axes = None

    def __init__(self, wcs):
        super(GWCSAdapter, self).__init__(wcs)

    def world_to_pixel(self, world_array):
        """
        Method for performing the world to pixel transformations.
        """
        return self.wcs.backward_transform(world_array)

    def pixel_to_world(self, pixel_array):
        """
        Method for performing the pixel to world transformations.
        """
        return self.wcs(pixel_array, with_bounding_box=False) * self.spectral_axis_unit

    @property
    def spectral_axis_unit(self):
        """
        Returns the unit of the spectral axis.
        """
        return self._wcs.output_frame.unit[0]

    @property
    def rest_frequency(self):
        """
        Returns the rest frequency defined in the WCS.
        """
        return

    @property
    def rest_wavelength(self):
        """
        Returns the rest wavelength defined in the WCS.
        """
        return
