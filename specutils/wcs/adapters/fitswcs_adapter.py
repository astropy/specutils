from astropy.wcs import (WCS, WCSSUB_CELESTIAL, WCSSUB_CUBEFACE,
                         WCSSUB_LATITUDE, WCSSUB_LONGITUDE, WCSSUB_SPECTRAL,
                         WCSSUB_STOKES)
from astropy.wcs import InvalidSubimageSpecificationError

# Use this once in specutils
from ...utils.wcs_utils import convert_spectral_axis, determine_ctype_from_vconv
from ..wcs_adapter import WCSAdapter, WCSAxes

__all__ = ['FITSWCSAdapter']


class FITSWCSAdapter(WCSAdapter):
    """
    Adapter class that adds support for FITSWCS objects.
    """
    wrapped_class = WCS
    axes = None

    def __init__(self, wcs):
        super(FITSWCSAdapter, self).__init__(wcs)

        # Store a reference to all axes information within the wcs object
        self.axes = WCSAxes(
            longitude=self.wcs.sub([WCSSUB_LONGITUDE]),
            latitude=self.wcs.sub([WCSSUB_LATITUDE]),
            cubeface=self.wcs.sub([WCSSUB_CUBEFACE]),
            spectral=self.wcs.sub([WCSSUB_SPECTRAL]),
            stokes=self.wcs.sub([WCSSUB_STOKES]),
            celestial=self.wcs.sub([WCSSUB_CELESTIAL])
        )

        # TODO: make this more efficient. Check to see whether the spectral
        # axis was actually parsed
        if self.axes.spectral.naxis == 0:
            try:
                idx = list(self.wcs.wcs.ctype).index('LINEAR')
                self.axes = self.axes._replace(spectral=self.wcs.sub([idx + 1]))
            except ValueError:
                raise InvalidSubimageSpecificationError(
                    "Cannot find a spectral axis in the provided WCS."
                    "Are your 'ctype's correct?")

    def world_to_pixel(self, world_array):
        """
        Method for performing the world to pixel transformations.
        """
        return self.axes.spectral.all_world2pix(world_array, 0)[0]

    def pixel_to_world(self, pixel_array):
        """
        Method for performing the pixel to world transformations.
        """
        return self.axes.spectral.all_pix2world(pixel_array, 0)[0]

    @property
    def spectral_axis_unit(self):
        """
        Returns the unit of the spectral axis.
        """
        return self._wcs.wcs.cunit[self._wcs.wcs.spec]

    @property
    def rest_frequency(self):
        """
        Returns the rest frequency defined in the WCS.
        """
        return self.wcs.wcs.restfrq

    @property
    def rest_wavelength(self):
        """
        Returns the rest wavelength defined in the WCS.
        """
        return self.wcs.wcs.restwav

    def bin_edges(self):
        # the WCS doesn't know about its own pixel array
        edge_indices = list(self.axes.spectral.pixel_indices - 0.5) + \
                       [self.axes.spectral.pixel_indices[-1] + 0.5]

        return self.pixel_to_world(edge_indices, 0)

    def with_spectral_unit(self, unit, rest_value, velocity_convention):
        # Shorter versions to keep lines under 80
        ctype_from_vconv = determine_ctype_from_vconv

        out_ctype = ctype_from_vconv(self._wcs.wcs.ctype[self._wcs.wcs.spec],
                                     unit,
                                     velocity_convention=velocity_convention)

        new_wcs = convert_spectral_axis(self._wcs, unit, out_ctype,
                                        rest_value=rest_value)

        new_wcs.wcs.set()

        return new_wcs
