import os

from ..interfaces.importer import import_modules

modules = import_modules(
    os.path.abspath(os.path.join(__file__, '..', 'loaders')))