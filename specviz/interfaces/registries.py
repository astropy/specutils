from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import yaml
import os
import importlib
import inspect
import logging
from functools import wraps

import astropy.io.registry as io_registry

from ..core.data import GenericSpectrum1D
from ..ui.widgets.plugin import Plugin


class Registry(object):
    """
    Maintains a set of referential objects.
    """
    def __init__(self):
        self._members = []

    @property
    def members(self):
        return self._members

    def add(self):
        raise NotImplementedError()

    def get(self):
        raise NotImplementedError()


class PluginRegistry(Registry):
    """Loads and stores custom plugins."""
    def __init__(self):
        super(PluginRegistry, self).__init__()

        cur_path = os.path.abspath(os.path.join(__file__, '..', '..',
                                                'plugins'))
        for mod in os.listdir(cur_path):
            mod = mod.split('.')[0]
            mod = importlib.import_module("specviz.plugins." + str(mod))
            cls_members = inspect.getmembers(
                mod, lambda member: inspect.isclass(member)
                                    and Plugin in member.__bases__)

            if len(cls_members) == 0:
                continue

            for _, cls_plugin in cls_members:
                self._members.append(cls_plugin())

    def add(self):
        pass


class LoaderRegistry(Registry):
    def __init__(self):
        super(LoaderRegistry, self).__init__()

    def __call__(self, label, identifier, priority=-1, **kwargs):
        def decorator(func):
            func.loader_wrapper = True

            @wraps(func)
            def wrapper(label, identifier, priority=-1, **kwargs):
                logging.info("Added {} to loader registry.".format(label))
                io_registry.register_reader(label, GenericSpectrum1D,
                                            func)
                io_registry.register_identifier(label, GenericSpectrum1D,
                                                identifier)

                return func

            return wrapper
        return decorator

    def _load_py(self):
        """
        Loads python files as custom loaders.
        """
        cur_path = os.path.abspath(os.path.join(__file__, '..', '..', 'io',
                                                'loaders'))

        for mod in os.listdir(cur_path):
            mod = mod.split('.')[0]
            mod = importlib.import_module("specviz.loaders." + str(mod))
            members = inspect.getmembers(
                mod, lambda member: inspect.isfunction(member))

            if len(members) == 0:
                continue

            for _, func in members:
                if func.loader_wrapper:
                    self._members.append(func)

    def _load_yaml(self):
        """
        Loads yaml files as custom loaders.
        """
        cur_path = os.path.join(os.path.dirname(__file__), 'yaml_loaders')
        usr_path = os.path.join(os.path.expanduser('~'), '.specviz')
        lines_path = os.path.join(os.path.dirname(__file__), '../data/linelists')

        # This order determines priority in case of duplicates; paths higher
        # in this list take precedence
        check_paths = [usr_path, cur_path, lines_path]

        if not os.path.exists(usr_path):
            os.mkdir(usr_path)

        for path in check_paths:
            for file_name in [x for x in os.listdir(path)
                              if x.endswith('yaml')]:
                f_path = os.path.join(path, file_name)
                custom_loader = yaml.load(open(f_path, 'r'))
                custom_loader.set_filter()

                self.add(custom_loader)

    def add(self, loader):
        if len([x for x in self._members if x.filter == loader.filter]) == 0:
            self._members.append(loader)

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
loader_registry = LoaderRegistry()
plugin_registry = PluginRegistry()
