import abc
import logging
from collections import namedtuple

WCSAxes = namedtuple('WCSAxes', ['longitude', 'latitude', 'cubeface',
                                 'spectral', 'stokes', 'celestial', ])


class WCSAdapterMeta(type):
    """
    Meta class to add any subclasses to a registry associating some WCS object
    with a specific ~`WCSAdapter` class that `specutils` can interface with.
    """
    def __init__(cls, name, bases, dct):
        if not hasattr(cls, 'registry'):
            cls.registry = {}
        else:
            wcs_class = dct.get('wrapped_class', None)

            if wcs_class is not None:
                logging.debug("Added {} to adapter registry.".format(wcs_class))
                cls.registry[wcs_class] = cls

        super(WCSAdapterMeta, cls).__init__(name, bases, dct)


class WCSAdapter(metaclass=type('WCSAdapterMetaProxy',
                                (abc.ABCMeta, WCSAdapterMeta), {})):
    """
    An adapter class that exposes the specific methods and attributes that
    `specutils` understands.
    """

    def __init__(self, wcs):
        """
        Initialize the adapter class.

        Parameters
        ----------
        wcs : object
            Any wcs object, the class of which will be stored in the registry
            in association with this adapter.
        """
        self._wcs = wcs

    @property
    def wcs(self):
        """
        Reference to the stored WCS object.

        Returns
        -------
        object
            The WCS object associated with this adapter class.
        """
        return self._wcs

    @abc.abstractmethod
    def axes(self):
        """
        The associated axes described in the WCS object.
        """
        pass

    @abc.abstractmethod
    def wrapped_class(self):
        """
        A specific reference to the WCS class to which this object should be
        associated.
        """
        pass

    @abc.abstractmethod
    def world_to_pixel(self, world_array):
        """
        Method for performing the world to pixel transformations.
        """
        pass

    @abc.abstractmethod
    def pixel_to_world(self, pixel_array):
        """
        Method for performing the pixel to world transformations.
        """
        pass

    @abc.abstractmethod
    def spectral_axis_unit(self):
        """
        Returns the unit of the spectral axis.
        """
        pass

    @abc.abstractmethod
    def rest_frequency(self):
        pass

    @abc.abstractmethod
    def rest_wavelength(self):
        pass

    @abc.abstractmethod
    def with_spectral_unit(self):
        pass

    @abc.abstractmethod
    def __getitem__(self, item):
        pass

    def __repr__(self):
        return "<Identity Transform WCS: pixel - {} transformation>".format(
            self.spectral_axis_unit)
