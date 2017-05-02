import six
import abc
from astropy.wcs import WCS


class WCSAdapterMeta(type):
    def __new__(cls, name, parents, dct):
        # create a class_id if it's not specified
        if 'class_id' not in dct:
            dct['class_id'] = name.lower()

        # open the specified file for writing
        if 'file' in dct:
            filename = dct['file']
            dct['file'] = open(filename, 'w')

        # we need to call type.__new__ to complete the initialization
        return super(WCSAdapterMeta, cls).__new__(cls, name, parents, dct)

    def __init__(cls, name, bases, dct):
        if not hasattr(cls, 'registry'):
            cls.registry = {}
        else:
            wcs_class = dct.get('wrapped_class', None)

            if wcs_class is None:
                raise NotImplemented

            # cls.registry.get(wcs_class, []).append(cls)
            cls.registry[wcs_class] = cls

        super(WCSAdapterMeta, cls).__init__(name, bases, dct)


class WCSAdapterMetaProxy(abc.ABCMeta, WCSAdapterMeta):
    pass


@six.add_metaclass(WCSAdapterMetaProxy)
class WCSAdapter:
    def __init__(self, wcs):
        self._wcs = wcs

    @abc.abstractproperty
    def axes(self):
        pass

    @abc.abstractmethod
    def world_to_pix(self):
        pass

    @abc.abstractmethod
    def pix_to_world(self):
        pass


class FITSWCSAdapter(WCSAdapter):
    wrapped_class = WCS

    axes = None

    def world_to_pix(self):
        pass

    def pix_to_world(self):
        pass