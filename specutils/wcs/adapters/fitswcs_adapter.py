from ..wcs_adapter import WCSAdapter, WCSAxes

from astropy.wcs import WCS
from astropy.wcs import WCSSUB_SPECTRAL, WCSSUB_LONGITUDE, WCSSUB_LATITUDE, \
    WCSSUB_CUBEFACE, WCSSUB_STOKES, WCSSUB_CELESTIAL

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

    def world_to_pixel(self, world_array):
        """
        Method for performing the world to pixel transformations.
        """
        return self.axes.spectral.all_world2pixel(world_array, 0)[0]

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
        return self.axes.spectral.wcs.cunit[0]

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
        edge_indices = list(self.axes.spectral.pixel_indices-0.5) + \
                       [self.axes.spectral.pixel_indices[-1]+0.5]

        return self.pix2world(edge_indices, 0)