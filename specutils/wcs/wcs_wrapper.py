from .wcs_adapter import WCSAdapter
from .adapters import *

import gwcs
from gwcs import coordinate_frames as cf
from astropy.modeling.tabular import Tabular1D
from astropy.units import quantity_input
import numpy as np

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
        if wcs is None:
            return

        adapter_class = WCSAdapter.registry.get(wcs.__class__, None)

        if adapter_class is not None:
            return adapter_class(wcs, *args, **kwargs)

        raise NotImplementedError("No such adapter for class "
                                  "{}".format(wcs.__class__))

    @staticmethod
    @quantity_input(array=['length', 'frequency', 'speed'])
    def from_array(array):
        """
        Create a new WCS from provided tabular data. This defaults to being
        a GWCS object.
        """
        coord_frame = cf.CoordinateFrame(naxes=1,
                                         axes_type=('SPECTRAL',),
                                         axes_order=(0,))
        spec_frame = cf.SpectralFrame(unit=array.unit, axes_order=(0,))
        transform = Tabular1D(np.arange(len(array)), array.value)

        tabular_gwcs = gwcs.wcs.WCS(
            [(coord_frame, transform), (spec_frame, None)])

        return WCSWrapper(wcs=tabular_gwcs)
