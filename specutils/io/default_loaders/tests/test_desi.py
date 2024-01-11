from astropy.io import fits
from astropy.utils.data import download_file
from astropy.config import set_temp_cache
import pytest

# Package
from ..desi import (coadd_loader, spectra_loader)
from ....spectra import SpectrumList


has_importlib = True
try:
    from importlib.resources import files
except ImportError:
    from pkg_resources import resource_filename as files
    has_importlib = False


@pytest.fixture(scope="function")
def local_filename(request):
    if has_importlib:
        return files('specutils.io.default_loaders.tests') / 't' / request.param
    else:
        return files('specutils.io.default_loaders.tests', f't/{request.param}')


@pytest.fixture(scope="function")
def remote_filename(request):
    if 'dark' in request.param:
        return f'https://data.desi.lbl.gov/public/edr/spectro/redux/fuji/healpix/sv3/dark/260/26065/{request.param}'
    else:
        return f'https://data.desi.lbl.gov/public/edr/spectro/redux/fuji/tiles/cumulative/169/20210419/{request.param}'


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
                                            'spectra-5-169-thru20210419.fits'], indirect=['local_filename'])
def test_SpectrumList_reader(local_filename):
    spectrum = SpectrumList.read(local_filename)  # noqa
    assert spectrum[0].meta['band'] == 'b'
    assert spectrum[1].meta['band'] == 'r'
    assert spectrum[2].meta['band'] == 'z'


@pytest.mark.parametrize('local_filename, fmt', [('coadd-sv3-dark-26065.fits', 'DESI coadd'),
                                                 ('spectra-sv3-dark-26065.fits', 'DESI spectra'),
                                                 ('coadd-5-169-thru20210419.fits', 'DESI coadd'),
                                                 ('spectra-5-169-thru20210419.fits', 'DESI spectra')], indirect=['local_filename'])
def test_SpectrumList_reader_fileobj(local_filename, fmt):
    with open(local_filename, 'rb') as fileobj:
        spectrum = SpectrumList.read(fileobj, format=fmt)  # noqa
    assert spectrum[0].meta['band'] == 'b'
    assert spectrum[1].meta['band'] == 'r'
    assert spectrum[2].meta['band'] == 'z'


@pytest.mark.parametrize('local_filename, fmt', [('coadd-sv3-dark-26065.fits', 'DESI coadd'),
                                                 ('spectra-sv3-dark-26065.fits', 'DESI spectra'),
                                                 ('coadd-5-169-thru20210419.fits', 'DESI coadd'),
                                                 ('spectra-5-169-thru20210419.fits', 'DESI spectra')], indirect=['local_filename'])
def test_SpectrumList_reader_hdulist(local_filename, fmt):
    with fits.open(local_filename, mode='readonly') as hdulist:
        spectrum = SpectrumList.read(hdulist, format=fmt)  # noqa
    assert spectrum[0].meta['band'] == 'b'
    assert spectrum[1].meta['band'] == 'r'
    assert spectrum[2].meta['band'] == 'z'


@pytest.mark.parametrize('remote_filename, loader', [pytest.param('coadd-sv3-dark-26065.fits', coadd_loader, marks=pytest.mark.remote_data),
                                                     pytest.param('spectra-sv3-dark-26065.fits', spectra_loader, marks=pytest.mark.remote_data),
                                                     pytest.param('coadd-5-169-thru20210419.fits', coadd_loader, marks=pytest.mark.remote_data),
                                                     pytest.param('spectra-5-169-thru20210419.fits', spectra_loader, marks=pytest.mark.remote_data)], indirect=['remote_filename'])
def test_loader_remote(tmp_path, remote_filename, loader):
    with set_temp_cache(path=str(tmp_path)):
        filename = download_file(remote_filename, cache=True)
        spectrum = loader(filename)  # noqa
        assert spectrum[0].meta['band'] == 'b'
        assert spectrum[1].meta['band'] == 'r'
        assert spectrum[2].meta['band'] == 'z'


@pytest.mark.parametrize('remote_filename, fmt', [pytest.param('coadd-sv3-dark-26065.fits', 'DESI coadd', marks=pytest.mark.remote_data),
                                                  pytest.param('spectra-sv3-dark-26065.fits', 'DESI spectra', marks=pytest.mark.remote_data),
                                                  pytest.param('coadd-5-169-thru20210419.fits', 'DESI coadd', marks=pytest.mark.remote_data),
                                                  pytest.param('spectra-5-169-thru20210419.fits', 'DESI spectra', marks=pytest.mark.remote_data)], indirect=['remote_filename'])
def test_SpectrumList_reader_remote(tmp_path, remote_filename, fmt):
    with set_temp_cache(path=str(tmp_path)):
        filename = download_file(remote_filename, cache=True)
        spectrum = SpectrumList.read(filename, format=fmt)  # noqa
        assert spectrum[0].meta['band'] == 'b'
        assert spectrum[1].meta['band'] == 'r'
        assert spectrum[2].meta['band'] == 'z'
