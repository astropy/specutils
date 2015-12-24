from __future__ import absolute_import, division, print_function

import logging
from astropy.io import registry as io_registry
from ..core.data import Data

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import StdDevUncertainty

import yaml
import os


class Registry(object):
    """
    Maintains a set of referential objects.
    """
    def __init__(self):
        self._members = []


class CustomLoaderRegistry(Registry):
    def __init__(self):
        super(CustomLoaderRegistry, self).__init__()

        cur_path = os.path.join(os.path.dirname(__file__), 'default_loaders')

        for file_name in os.listdir(cur_path):
            f_path = os.path.join(cur_path, file_name)
            custom_loader = yaml.load(open(f_path, 'r'))
            custom_loader.set_filter()

            self._members.append(custom_loader)

    def get(self, filter):
        return [x for x in self._members if x.filter == filter][0]

    @property
    def filters(self):
        return [x.filter for x in self._members]


class YAMLLoader(yaml.YAMLObject):
    yaml_tag = u'!CustomLoader'

    def __init__(self, extension, name, data, uncertainty, mask, meta, wcs):
        self.name = name
        self.extension = extension
        self.data = data
        self.uncertainty = uncertainty
        self.mask = mask
        self.meta = meta
        self.wcs = wcs
        self.filter = None

    def set_filter(self):
        if isinstance(self.extension, list):
            filter_string = ' '.join(['*.{}'.format(x)
                                       for x in self.extension])
            self.filter = "{} ({})".format(self.name, filter_string)
        else:
            self.filter = "{} (*.{})".format(self.name, self.extension)


# Create loader registry instance
loader_registry = CustomLoaderRegistry()


def fits_reader(filename, filter, **kwargs):
    """
    This generic function will query the loader factory, which has already
    loaded the yaml configuration files, in an attempt to parse the
    associated fits file.
    """
    hdulist = fits.open(filename, **kwargs)
    ref = loader_registry.get(filter)

    wcs = WCS(hdulist[ref.wcs['hdu']].header)
    data = hdulist[ref.data['hdu']].data

    if ref.data.get('col') is not None:
        data = data[data.columns[ref.data['col']].name]

    uncertainty = hdulist[ref.uncertainty['hdu']].data
    uncertainty_type = ref.uncertainty.get('type') or 'var'

    uncertainty = StdDevUncertainty(uncertainty)

    mask = np.zeros(shape=data.shape)

    try:
        if ref.mask.get('hdu') is not None:
            mask = hdulist[ref.mask['hdu']].data
    except IndexError:
        logging.warning("Incorrect mask HDU; ignoring.")

    return Data(data=data, uncertainty=uncertainty, mask=mask, wcs=wcs)


def fits_identify(origin, *args, **kwargs):
    return isinstance(args[0], basestring) and \
           args[0].lower().split('.')[-1] in ['fits', 'fit']


# Add IO reader/identifier to io registry
io_registry.register_reader('fits', Data, fits_reader)
io_registry.register_identifier('fits', Data, fits_identify)