import os

from asdf.util import filepath_to_url
from asdf.extension import AsdfExtension

from astropy.io.misc.asdf.extension import ASTROPY_SCHEMA_URI_BASE

from .tags.spectra import *
from .types import _specutils_types


SCHEMA_PATH = os.path.abspath(
    os.path.join(os.path.dirname(__file__), 'schemas'))
SPECUTILS_URL_MAPPING = [
    (ASTROPY_SCHEMA_URI_BASE,
     filepath_to_url(
         os.path.join(SCHEMA_PATH, 'astropy.org')) +
         '/{url_suffix}.yaml')]


class SpecutilsExtension(AsdfExtension):
    @property
    def types(self):
        return _specutils_types

    @property
    def tag_mapping(self):
        return [('tag:astropy.org:specutils',
                ASTROPY_SCHEMA_URI_BASE + 'specutils{tag_suffix}')]

    @property
    def url_mapping(self):
        return SPECUTILS_URL_MAPPING
