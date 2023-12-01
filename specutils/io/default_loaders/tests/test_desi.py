from astropy.utils.data import download_file
from astropy.config import set_temp_cache
import pytest

# Package
from ..desi import (coadd_loader, spectra_loader)


@pytest.mark.remote_data
def test_coadd_loader(tmp_path):
    coadd = 'https://data.desi.lbl.gov/public/edr/spectro/redux/fuji/healpix/sv3/dark/260/26065/coadd-sv3-dark-26065.fits'
    with set_temp_cache(path=str(tmp_path)):
        filename = download_file(coadd, cache=True)
        spectrum = coadd_loader(filename)  # noqa
        assert spectrum[0].meta['band'] == 'b'
        assert spectrum[1].meta['band'] == 'r'
        assert spectrum[2].meta['band'] == 'z'


@pytest.mark.remote_data
def test_spectra_loader(tmp_path):
    spectra = 'https://data.desi.lbl.gov/public/edr/spectro/redux/fuji/healpix/sv3/dark/260/26065/spectra-sv3-dark-26065.fits'
    with set_temp_cache(path=str(tmp_path)):
        filename = download_file(spectra, cache=True)
        spectrum = spectra_loader(filename)  # noqa
        assert spectrum[0].meta['band'] == 'b'
        assert spectrum[1].meta['band'] == 'r'
        assert spectrum[2].meta['band'] == 'z'
