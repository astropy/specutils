from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import yaml
import os
from astropy.io import registry as io_registry


class Registry(object):
    """
    Maintains a set of referential objects.
    """
    def __init__(self):
        self._members = []


class CustomLoaderRegistry(Registry):
    def __init__(self):
        super(CustomLoaderRegistry, self).__init__()

        check_paths = []

        cur_path = os.path.join(os.path.dirname(__file__), 'default_loaders')
        usr_path = os.path.join(os.path.expanduser('~'), '.pyfocal')

        if not os.path.exists(usr_path):
            os.mkdir(usr_path)

        for path in check_paths:
            for file_name in os.listdir(path):
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

    def __init__(self, extension, name, data, dispersion, uncertainty, mask,
                 wcs, meta):
        self.name = name
        self.extension = extension
        self.data = data
        self.dispersion = dispersion
        self.uncertainty = uncertainty
        self.mask = mask
        self.wcs = wcs
        self.meta = meta or {}
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

# Import loaders
from .loaders import register_loaders
register_loaders()
