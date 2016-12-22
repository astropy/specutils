"""
Registry library
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import yaml
import os, sys
import importlib
import inspect
import logging

import astropy.io.registry as io_registry

#-- local
from ..core.data import Spectrum1DRef
from ..io.yaml_loader import FitsYamlRegister, AsciiYamlRegister
from ..io.loaders import *

__all__ = ['Registry',
           'PluginRegistry',
           'LoaderRegistry']


class Registry(object):
    """
    Maintains a set of referential objects.

    Attributes
    ----------
    members: list
        The list of members belonging to this registry.
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
        usr_path = os.path.join(os.path.expanduser('~'), '.specviz')

        # This order determines priority in case of duplicates; paths higher
        # in this list take precedence
        check_paths = [usr_path, cur_path]

        if not os.path.exists(usr_path):
            os.mkdir(usr_path)

        for path in check_paths:
            for mod in [x for x in os.listdir(path) if x.endswith('.py')]:
                mod = mod.split('.')[0]

                if path == cur_path:
                    mod = "specviz.plugins.{}".format(mod)
                else:
                    sys.path.insert(0, path)

                mod = importlib.import_module(mod)

                cls_members = inspect.getmembers(
                    mod, lambda member: inspect.isclass(member)
                                        and 'Plugin' in [x.__name__
                                                         for x in
                                                         member.__bases__])

                for cls_name, cls_plugin in cls_members:
                    self._members.append(cls_plugin())

                sys.path.pop(0)


class LoaderRegistry(Registry):
    """Loads and stores the IO data loaders
    """
    def __init__(self):
        super(LoaderRegistry, self).__init__()

        self._load_py()
        self._load_yaml()

    def _load_py(self):
        """ Loads built-in and custom python loaders

        Loaders from the io.loaders module will be included from the
        import statement.
        Python modules (.py ending) found in the following locations will be
        auto-loaded into the registry for data loading.

        1.  .specviz folder in the user's HOME directory

        """

        usr_path = os.path.join(os.path.expanduser('~'), '.specviz')

        # This order determines priority in case of duplicates; paths higher
        # in this list take precedence
        #
        # Leaving in list format incase other locations want to be added
        # in the future
        check_paths = [usr_path]

        if not os.path.exists(usr_path):
            os.mkdir(usr_path)

        for path in check_paths:
            for mod in [x for x in os.listdir(path) if x.endswith('.py')]:
                mod = mod.split('.')[0]
                sys.path.insert(0, path)
                mod = importlib.import_module(mod)

                # for _, func in members:
                #     if hasattr(func, 'loader_wrapper') and func.loader_wrapper:
                #         self._members.append(func)

                sys.path.pop(0)

    def _load_yaml(self):
        """ Loads yaml files as custom loaders.

        YAML files found in the following three locations will be auto-loaded
        into the registry for data loading.

        1.  .specviz folder in the user's HOME directory
        2.  the current working directory
        3.  the linelists directory delivered with this package.

        The io_registry will be updated with the YAML schematics for each of the
        different filetypes.  Errors in loading the registry will write an error
        to the log.

        """
        cur_path = os.path.join(os.path.dirname(__file__), '..', 'io',
                                'yaml_loaders')
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

                # Figure out which of the two generic loaders to associate
                # this yaml file with
                if any(ext in custom_loader.extension for ext in ['fits']):
                    loader = FitsYamlRegister(custom_loader)
                elif any(ext in custom_loader.extension
                         for ext in ['txt', 'data']):
                    loader = AsciiYamlRegister(custom_loader)

                try:
                    io_registry.register_reader(custom_loader.name,
                                                Spectrum1DRef,
                                                loader.reader)
                    io_registry.register_identifier(custom_loader.name,
                                                    Spectrum1DRef,
                                                    loader.identify)
                except io_registry.IORegistryError as e:
                    logging.error(e)


class YAMLLoader(yaml.YAMLObject):
    """ Helper to load YAML files
    """
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
