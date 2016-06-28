from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import importlib, inspect
import os
import glob, imp
from os.path import join, basename, splitext


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

