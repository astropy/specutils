# This file is used to configure the behavior of pytest when using the Astropy
# test infrastructure. It needs to live inside the package in order for it to
# get picked up when running the tests inside an interpreter using
# packagename.test

try:
    from pytest_astropy_header.display import PYTEST_HEADER_MODULES, TESTED_VERSIONS
    ASTROPY_HEADER = True
except ImportError:
    ASTROPY_HEADER = False


def pytest_configure(config):

    if ASTROPY_HEADER:

        config.option.astropy_header = True

        # Customize the following lines to add/remove entries from the list of
        # packages for which version numbers are displayed when running the tests.
        PYTEST_HEADER_MODULES.pop('Pandas', None)
        PYTEST_HEADER_MODULES['gwcs'] = 'gwcs'
        del PYTEST_HEADER_MODULES['h5py']
        del PYTEST_HEADER_MODULES['Matplotlib']
        # Use ASDF schema tester plugin if ASDF is installed
        from importlib.util import find_spec
        if find_spec('asdf') is not None:
            PYTEST_HEADER_MODULES['Asdf'] = 'asdf'

        from specutils import __version__
        TESTED_VERSIONS['specutils'] = __version__
