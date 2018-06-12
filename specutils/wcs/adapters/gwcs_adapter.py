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

    def __init__(self, wcs, unit=None):
        super(GWCSAdapter, self).__init__(wcs)

        self._rest_frequency = 0
        self._rest_wavelength = 0

        # TODO: Currently, unsure of how to create a copy of an arbitrary gwcs
        # object. For now, store the desired spectral axis unit information in
        # the adapter object instead of creating a new gwcs object like we do
        # for fitswcs.
        self._wcs_unit = self.wcs.output_frame.unit[0]
        self._unit = self._wcs_unit

        if unit is not None and unit.is_equivalent(self._unit,
                                                   equivalencies=u.spectral()):
            self._unit = unit

    def world_to_pixel(self, world_array):
        """
        Method for performing the world to pixel transformations.
        """
        return u.Quantity(self.wcs.invert(world_array),
                          self._wcs_unit).to(
            self._unit, equivalencies=u.spectral()).value

    def pixel_to_world(self, pixel_array):
        """
        Method for performing the pixel to world transformations.
        """
        return u.Quantity(self.wcs(pixel_array, with_bounding_box=False),
                          self._wcs_unit).to(
            self._unit, equivalencies=u.spectral()).value

    @property
    def spectral_axis_unit(self):
        """
        Returns the unit of the spectral axis.
        """
        return self._unit

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

    def with_spectral_unit(self, unit, rest_value=None, velocity_convention=None):
        """
        """
        if isinstance(unit, u.Unit) and unit.is_equivalent(
                self._wcs_unit, equivalencies=u.spectral()):
            return self.__class__(self.wcs, unit=unit)

        logging.error("WCS units incompatible: {} and {}.".format(
            unit, self._wcs_unit))
