from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import importlib, inspect
import os
import glob, imp
from os.path import join, basename, splitext


def import_plugins(path, filt_func):
    """
    Search a directory for importable `Plugin` classes.

    Parameters
    ----------
    path : str
        Path to directory containing python files with `Plugin` subclasses.
    filt_func : lambda
        A lambda function that filters the inspected python files.

    Returns
    -------
    instance_plugins : list
        A list of instantiated plugins.
    """
    instance_plugins = []

    for mod in os.listdir(path):
        mod = mod.split('.')[0]
        mod = importlib.import_module("specviz.plugins." + str(mod))
        cls_members = inspect.getmembers(mod, filt_func)

        if len(cls_members) == 0:
            continue

        for _, cls_plugin in cls_members:
            instance_plugins.append(cls_plugin())

    return instance_plugins


def import_modules(directory):
    """
    Import modules outside of this package into the current namespace.

    Parameters
    ----------
    directory : str
        Path to directory containing modules to be imported.

    Returns
    -------
    modules : dict
        A dictionary of modules where the key is the name and the value is
        the module object.
    """
    modules = {}

    for path in glob.glob(join(directory, '[!_]*.py')):
        name, ext = splitext(basename(path))
        modules[name] = imp.load_source(name, path)

    return modules

