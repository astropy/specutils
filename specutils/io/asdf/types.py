from asdf.types import CustomType, ExtensionTypeMeta


_specutils_types = set()


class SpecutilsTypeMeta(ExtensionTypeMeta):
    """
    Keeps track of `SpecutilsType` subclasses that are created so that they can
    be stored automatically by specutils extensions for ASDF.
    """
    def __new__(mcls, name, bases, attrs):
        cls = super().__new__(mcls, name, bases, attrs)
        # Classes using this metaclass are automatically added to the list of
        # specutils extensions
        _specutils_types.add(cls)
        return cls


class SpecutilsType(CustomType, metaclass=SpecutilsTypeMeta):
    """
    Parent class of all specutils tag implementations used by ASDF
    """
    organization = 'astropy.org'
    standard = 'specutils'
