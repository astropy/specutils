import six
import abc
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
                cls.registry[wcs_class] = cls

        super(WCSAdapterMeta, cls).__init__(name, bases, dct)


@six.add_metaclass(type('WCSAdapterMetaProxy',
                        (abc.ABCMeta, WCSAdapterMeta), {}))
class WCSAdapter:
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

    @abc.abstractproperty
    def axes(self):
        """
        The associated axes described in the WCS object.
        """
        pass

    @abc.abstractproperty
    def wrapped_class(self):
        """
        A specific reference to the WCS class to which this object should be
        associated.
        """
        pass

    @abc.abstractmethod
    def world_to_pix(self):
        """
        Method for performing the world to pixel transformations.
        """
        pass

    @abc.abstractmethod
    def pix_to_world(self):
        """
        Method for performing the pixel to world transformations.
        """
        pass