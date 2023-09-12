try:
    from pytest_astropy_header.display import PYTEST_HEADER_MODULES, TESTED_VERSIONS
    ASTROPY_HEADER = True
except ImportError:
    ASTROPY_HEADER = False


# Repeat this from specutils/conftest.py so tox picks it up.
def pytest_configure(config):

    if ASTROPY_HEADER:
        config.option.astropy_header = True

        # Customize the following lines to add/remove entries from the list of
        # packages for which version numbers are displayed when running the tests.
        PYTEST_HEADER_MODULES.pop('Pandas', None)
        PYTEST_HEADER_MODULES.pop('h5py', None)
        PYTEST_HEADER_MODULES['astropy'] = 'astropy'
        PYTEST_HEADER_MODULES['gwcs'] = 'gwcs'
        PYTEST_HEADER_MODULES['asdf'] = 'asdf'
        PYTEST_HEADER_MODULES['asdf-astropy'] = 'asdf_astropy'
        PYTEST_HEADER_MODULES['stdatamodels'] = 'stdatamodels'
        PYTEST_HEADER_MODULES['ndcube'] = 'ndcube'
        PYTEST_HEADER_MODULES['spectral-cube'] = 'spectral_cube'

        from specutils import __version__
        TESTED_VERSIONS['specutils'] = __version__
