import logging
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
        in the dispersion axis, this is just need to shift to the correct
        starting element.

        When WCS gets the abillity to do slicing then we might be able to
        remove this code.
        """
        # Create shift of x-axis
        if isinstance(item, int):
            shift = item
        elif isinstance(item, slice):
            shift = item.start
        else:
            raise TypeError('Unknown index type {}, must be int or slice.'.format(item))

        shifter = models.Shift(shift)

        # Get the current forward transform
        forward = self._wcs.forward_transform

        # Set the new transform
        self._wcs.set_transform(self._wcs.input_frame, self._wcs.output_frame, shifter|forward)
        """
        A container for `numpy.ndarray`-based datasets, using the
        `~astropy.nddata.NDDataBase` interface.

        The key distinction from raw `numpy.ndarray` is the presence of
        additional metadata such as uncertainty, mask, unit, a coordinate system
        and/or a dictionary containing further meta information. This class *only*
        provides a container for *storing* such datasets. For further functionality
        take a look at the ``See also`` section.

        Parameters
        -----------
        data : `numpy.ndarray`-like or `NDData`-like
            The dataset.

        mask : any type, optional
            Mask for the dataset. Masks should follow the ``numpy`` convention that
            **valid** data points are marked by ``False`` and **invalid** ones with
            ``True``.
            Defaults to ``None``.

        copy : `bool`, optional
            Indicates whether to save the arguments as copy. ``True`` copies
            every attribute before saving it while ``False`` tries to save every
            parameter as reference.
            Note however that it is not always possible to save the input as
            reference.
            Default is ``False``.

            .. versionadded:: 1.2

        Raises
        ------
        TypeError
            In case ``data`` or ``meta`` don't meet the restrictions.

        Notes
        -----
        Each attribute can be accessed through the homonymous instance attribute:
        ``data`` in a `NDData` object can be accessed through the `data`
        attribute::

            >>> from astropy.nddata import NDData
            >>> nd = NDData([1,2,3])
            >>> nd.data
            array([1, 2, 3])

        Given a conflicting implicit and an explicit parameter during
        initialization, for example the ``data`` is a `~astropy.units.Quantity` and
        the unit parameter is not ``None``, then the implicit parameter is replaced
        (without conversion) by the explicit one and a warning is issued::

            >>> import numpy as np
            >>> import astropy.units as u
            >>> q = np.array([1,2,3,4]) * u.m
            >>> nd2 = NDData(q, unit=u.cm)
            INFO: overwriting Quantity's current unit with specified unit. [astropy.nddata.nddata]
            >>> nd2.data  # doctest: +FLOAT_CMP
            array([1., 2., 3., 4.])
            >>> nd2.unit
            Unit("cm")

        See also
        --------
        NDDataArray
        """


        return self

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
        if isinstance(unit, u.UnitBase) and unit.is_equivalent(
                self._wcs_unit, equivalencies=u.spectral()):
            return self.__class__(self.wcs, unit=unit)

        logging.error("WCS units incompatible: {} and {}.".format(
            unit, self._wcs_unit))
