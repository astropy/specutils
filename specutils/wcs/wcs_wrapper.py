from .wcs_adapter import WCSAdapter
from .adapters import *


__all__ = ['WCSWrapper']


class WCSWrapper:
    """
    Factory class that returns an instance of a registered WCS adapter
    dependent upon the class of the `wcs` object passed in.
    """
    def __new__(cls, wcs, *args, **kwargs):
        """
        Control the creation of a new instance of a ~`WCSAdapter` subclass.

        Parameters
        ----------
        wcs : object
            A WCS object with a corresponding adapter in the WCS adapter
            registry.

        Returns
        -------
        object
            An instance of a `WCSAdapter` subclass that understands the WCS
            object used for initialization.
        """
        adapter_class = WCSAdapter.registry.get(wcs.__class__, None)

        if adapter_class is not None:
            return adapter_class(wcs, *args, **kwargs)

        raise NotImplementedError("No such adapater for class "
                                  "{}".format(wcs.__class__))