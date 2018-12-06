import logging
import copy

from gwcs.wcs import WCS
import astropy.units as u
from astropy.modeling import models

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

    def __getitem__(self, item):
        """
        This is a bit of a hack in order to fix the slicing of the WCS
        in the spectral dispersion direction.  The NDData slices properly
        but the spectral dispersion result was not.

        There is code slightly downstream that sets the *number* of entries
        in the dispersion axis, this is just needed to shift to the correct
        starting element.

        When WCS gets the ability to do slicing then we might be able to
        remove this code.
        """

        # Create shift of x-axis
        if isinstance(item, int):
            shift = item
        elif isinstance(item, slice):
            shift = item.start
        else:
            raise TypeError('Unknown index type {}, must be int or slice.'.format(item))

        # Create copy as we need to modify this and return it.
        newwcs = copy.deepcopy(self)

        if shift == 0:
            return newwcs

        shifter = models.Shift(shift)

        # Get the current forward transform
        forward = newwcs._wcs.forward_transform

        # Set the new transform
        newwcs.wcs.set_transform(newwcs._wcs.input_frame, newwcs._wcs.output_frame, shifter|forward)

        return newwcs

    def world_to_pixel(self, world_array):
        """
        Method for performing the world to pixel transformations.
        """
        # Convert world array to the wcs unit for the conversion
        with u.set_enabled_equivalencies(u.spectral()):
            world_array = u.Quantity(world_array, unit=self._wcs_unit)

        return self.wcs.invert(world_array.value)

    def pixel_to_world(self, pixel_array):
        """
        Method for performing the pixel to world transformations.
        """
        return u.Quantity(self.wcs(pixel_array, with_bounding_box=False),
                          self._wcs_unit).to(
            self.spectral_axis_unit, equivalencies=u.spectral())

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
        logging.debug("GWCS does not store rest frequency information. "
                      "Please define the rest value explicitly in the "
                      "`Spectrum1D` object.")

        return self._rest_frequency

    @property
    def rest_wavelength(self):
        """
        Returns the rest wavelength defined in the WCS.
        """
        logging.debug("GWCS does not store rest wavelength information. "
                      "Please define the rest value explicitly in the "
                      "`Spectrum1D` object.")

        return self._rest_wavelength

    def with_spectral_unit(self, unit, rest_value=None, velocity_convention=None):
        """
        """
        if isinstance(unit, u.UnitBase) and unit.is_equivalent(
                self._wcs_unit, equivalencies=u.spectral()):
            return self.__class__(self.wcs, unit=unit)

        logging.error("WCS units incompatible: {} and {}.".format(
            unit, self._wcs_unit))
