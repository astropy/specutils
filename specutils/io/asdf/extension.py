"""
Defines extension that is used by ASDF for recognizing specutils types
"""
import os
import urllib
import warnings

from asdf.util import filepath_to_url
from asdf.extension import AsdfExtension
from asdf.exceptions import AsdfDeprecationWarning

from .tags.spectra import *  # noqa
from .types import _specutils_types


with warnings.catch_warnings():
    warnings.filterwarnings(
        "ignore",
        category=AsdfDeprecationWarning,
        message=r".*from astropy.io.misc.asdf.* subclasses the deprecated CustomType .*",
    )
    from astropy.io.misc.asdf.extension import ASTROPY_SCHEMA_URI_BASE


SCHEMA_PATH = os.path.abspath(
    os.path.join(os.path.dirname(__file__), 'schemas'))
SPECUTILS_URL_MAPPING = [
    (urllib.parse.urljoin(ASTROPY_SCHEMA_URI_BASE, 'specutils/'),
     filepath_to_url(
         os.path.join(SCHEMA_PATH, 'astropy.org', 'specutils')) +
     '/{url_suffix}.yaml')]


class SpecutilsExtension(AsdfExtension):
    """
    Defines specutils types and schema locations to be used by ASDF
    """
    @property
    def types(self):
        """
        Collection of tag types that are used by ASDF for serialization
        """
        return _specutils_types

    @property
    def tag_mapping(self):
        """
        Defines mapping of specutils tag URIs to URLs
        """
        return [('tag:astropy.org:specutils',
                ASTROPY_SCHEMA_URI_BASE + 'specutils{tag_suffix}')]

    @property
    def url_mapping(self):
        """
        Defines mapping of specutils schema URLs into real locations on disk
        """
        return SPECUTILS_URL_MAPPING
