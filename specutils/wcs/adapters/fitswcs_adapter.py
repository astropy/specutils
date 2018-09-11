import astropy.units as u
from astropy.wcs import (WCS, WCSSUB_CELESTIAL, WCSSUB_CUBEFACE,
                         WCSSUB_LATITUDE, WCSSUB_LONGITUDE, WCSSUB_SPECTRAL,
                         WCSSUB_STOKES, InvalidSubimageSpecificationError)

# Use this once in specutils
from ...utils.wcs_utils import (convert_spectral_axis,
                                determine_ctype_from_vconv)
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
        self._spec_axis = None

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
            self.axes = self.axes._replace(spectral=self.wcs.sub([self.spec_axis + 1]))

    def __getitem__(self, item):
        """Pass slicing information to the internal `FITSWCS` object."""
        return self.wcs[item]

    def world_to_pixel(self, world_array):
        """
        Method for performing the world to pixel transformations.
        """
        return self.axes.spectral.all_world2pix(world_array, 0)[0]

    def pixel_to_world(self, pixel_array):
        """
        Method for performing the pixel to world transformations.
        """
        return u.Quantity(self.axes.spectral.all_pix2world(pixel_array, 0)[0],
                          self.spectral_axis_unit)

    @property
    def spec_axis(self):
        """
        Try and parse the spectral axis of the fits wcs object.
        """
        self._spec_axis = self.wcs.wcs.spec

        if self._spec_axis < 0:
            try:
                idx = list(self.wcs.wcs.ctype).index('LINEAR')
            except ValueError:
                raise InvalidSubimageSpecificationError(
                    "Cannot find a spectral axis in the provided WCS."
                    "Are your 'ctype's correct?")

            if self._wcs.wcs.spec < 0:
                self._spec_axis = idx

        return self._spec_axis

    @property
    def spectral_axis_unit(self):
        """
        Returns the unit of the spectral axis.
        """
        return self._wcs.wcs.cunit[self.spec_axis]

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

        out_ctype = ctype_from_vconv(self._wcs.wcs.ctype[self.spec_axis],
                                     unit,
                                     velocity_convention=velocity_convention)

        new_wcs = convert_spectral_axis(self._wcs, unit, out_ctype,
                                        rest_value=rest_value)

        new_wcs.wcs.set()

        return new_wcs
