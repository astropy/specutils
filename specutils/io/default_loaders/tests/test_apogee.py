# Third-party
from astropy.utils.data import download_file
from astropy.config import set_temp_cache
import pytest

# Package
from specutils.io.default_loaders.apogee import (apStar_loader, apVisit_loader,
                                                 aspcapStar_loader)


@pytest.mark.remote_data
def test_apStar_loader(tmp_path):
    apstar_url = ("https://data.sdss.org/sas/dr16/apogee/spectro/redux/r12/"
                  "stars/apo25m/N7789/apStar-r12-2M00005414+5522241.fits")
    with set_temp_cache(path=str(tmp_path)):
        filename = download_file(apstar_url, cache=True)
        spectrum = apStar_loader(filename)  # noqa


@pytest.mark.remote_data
def test_apVisit_loader(tmp_path):
    apvisit_url = ("https://data.sdss.org/sas/dr16/apogee/spectro/redux/r12/"
                   "visit/apo25m/N7789/5094/55874/"
                   "apVisit-r12-5094-55874-123.fits")
    with set_temp_cache(path=str(tmp_path)):
        filename = download_file(apvisit_url, cache=True)
        spectrum = apVisit_loader(filename)  # noqa


@pytest.mark.remote_data
def test_aspcapStar_loader(tmp_path):
    aspcap_url = ("https://data.sdss.org/sas/dr16/apogee/spectro/aspcap/r12/"
                  "l33/apo25m/N7789/aspcapStar-r12-2M00005414+5522241.fits")
    with set_temp_cache(path=str(tmp_path)):
        filename = download_file(aspcap_url, cache=True)
        spectrum = aspcapStar_loader(filename)  # noqa
