has_importlib = True
try:
    from importlib.resources import files
except ImportError:
    from pkg_resources import resource_filename as files
    has_importlib = False
from astropy.utils.data import download_file
from astropy.config import set_temp_cache
import pytest

# Package
from ..desi import (coadd_loader, spectra_loader)
from ....spectra import SpectrumList


@pytest.fixture(scope="function")
def local_filename(request):
    if has_importlib:
        return files('specutils.io.default_loaders.tests') / 't' / request.param
    else:
        return files('specutils.io.default_loaders.tests', f't/{request.param}')


@pytest.mark.parametrize('local_filename, loader', [('coadd-sv3-dark-26065.fits', coadd_loader),
                                                    ('spectra-sv3-dark-26065.fits', spectra_loader),
                                                    ('coadd-5-169-thru20210419.fits', coadd_loader),
                                                    ('spectra-5-169-thru20210419.fits', spectra_loader)], indirect=['local_filename'])
def test_loader(local_filename, loader):
    spectrum = loader(local_filename)  # noqa
    assert spectrum[0].meta['band'] == 'b'
    assert spectrum[1].meta['band'] == 'r'
    assert spectrum[2].meta['band'] == 'z'


@pytest.mark.parametrize('local_filename', ['coadd-sv3-dark-26065.fits',
                                            'spectra-sv3-dark-26065.fits',
                                            'coadd-5-169-thru20210419.fits',
                                            'spectra-5-169-thru20210419.fits'], indirect=True)
def test_SpectrumList_reader(local_filename):
    spectrum = SpectrumList.read(local_filename)  # noqa
    assert spectrum[0].meta['band'] == 'b'
    assert spectrum[1].meta['band'] == 'r'
    assert spectrum[2].meta['band'] == 'z'


@pytest.mark.remote_data
def test_coadd_loader_remote(tmp_path):
    coadd = 'https://data.desi.lbl.gov/public/edr/spectro/redux/fuji/healpix/sv3/dark/260/26065/coadd-sv3-dark-26065.fits'
    with set_temp_cache(path=str(tmp_path)):
        filename = download_file(coadd, cache=True)
        spectrum = coadd_loader(filename)  # noqa
        assert spectrum[0].meta['band'] == 'b'
        assert spectrum[1].meta['band'] == 'r'
        assert spectrum[2].meta['band'] == 'z'


@pytest.mark.remote_data
def test_spectra_loader_remote(tmp_path):
    spectra = 'https://data.desi.lbl.gov/public/edr/spectro/redux/fuji/healpix/sv3/dark/260/26065/spectra-sv3-dark-26065.fits'
    with set_temp_cache(path=str(tmp_path)):
        filename = download_file(spectra, cache=True)
        spectrum = spectra_loader(filename)  # noqa
        assert spectrum[0].meta['band'] == 'b'
        assert spectrum[1].meta['band'] == 'r'
        assert spectrum[2].meta['band'] == 'z'


@pytest.mark.remote_data
def test_SpectrumList_coadd_reader_remote(tmp_path):
    coadd = 'https://data.desi.lbl.gov/public/edr/spectro/redux/fuji/healpix/sv3/dark/260/26065/coadd-sv3-dark-26065.fits'
    with set_temp_cache(path=str(tmp_path)):
        filename = download_file(coadd, cache=True)
        spectrum = SpectrumList.read(filename)  # noqa
        assert spectrum[0].meta['band'] == 'b'
        assert spectrum[1].meta['band'] == 'r'
        assert spectrum[2].meta['band'] == 'z'


@pytest.mark.remote_data
def test_SpectrumList_spectra_reader_remote(tmp_path):
    spectra = 'https://data.desi.lbl.gov/public/edr/spectro/redux/fuji/healpix/sv3/dark/260/26065/spectra-sv3-dark-26065.fits'
    with set_temp_cache(path=str(tmp_path)):
        filename = download_file(spectra, cache=True)
        spectrum = SpectrumList.read(filename)  # noqa
        assert spectrum[0].meta['band'] == 'b'
        assert spectrum[1].meta['band'] == 'r'
        assert spectrum[2].meta['band'] == 'z'
