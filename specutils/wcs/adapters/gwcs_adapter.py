import logging
from gwcs.wcs import WCS
import astropy.units as u

from ..wcs_adapter import WCSAdapter

__all__ = ['GWCSAdapter']


class GWCSAdapter(WCSAdapter):
    """
    Adapter class that adds support for GWCS objects.
    """
    wrapped_class = WCS
    axes = None

    def __init__(self, wcs):
        super(GWCSAdapter, self).__init__(wcs)

        self._rest_frequency = 0
        self._rest_wavelength = 0

    def world_to_pixel(self, world_array):
        """
        Method for performing the world to pixel transformations.
        """
        return self.wcs.invert(world_array)

    def pixel_to_world(self, pixel_array):
        """
        Method for performing the pixel to world transformations.
        """
        return self.wcs(pixel_array, with_bounding_box=False)

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
        logging.warning("GWCS does not store rest frequency information. "
                        "Please define the rest value explicitly in the "
                        "`Spectrum1D` object.")

        return self._rest_frequency

    @property
    def rest_wavelength(self):
        """
        Returns the rest wavelength defined in the WCS.
        """
        logging.warning("GWCS does not store rest wavelength information. "
                        "Please define the rest value explicitly in the "
                        "`Spectrum1D` object.")

        return self._rest_wavelength

    @property
    def unit(self):
        return self._wcs.unit

    def with_new_unit(self, unit, rest_value, velocity_convention):
        """
        """
        if isinstance(unit, u.Unit):
            new_output_frame = self._wcs.output_frame
            new_output_frame.unit = unit

            wcs = WCS(forward_transform=self._wcs.forward_transform,
                      input_frame=self._wcs.input_frame,
                      output_frame=new_output_frame)

            return wcs
