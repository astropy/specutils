# Third-party
import pytest
from astropy.config import set_temp_cache
from astropy.utils.data import download_file

# Package
from specutils.io.default_loaders.apogee import apStar_loader, apVisit_loader, aspcapStar_loader


@pytest.mark.remote_data
@pytest.mark.parametrize(
    "url",
    [
        "https://data.sdss.org/sas/dr16/apogee/spectro/redux/r12/stars/apo25m/N7789/apStar-r12-2M00005414+5522241.fits",
        "https://data.sdss.org/sas/dr17/apogee/spectro/redux/dr17/stars/lco25m/293-37-C/asStar-dr17-2M03174068-7748459.fits",
        "https://data.sdss.org/sas/dr13/apogee/spectro/redux/r6/stars/apo25m/4291/apStar-r6-2M07440053+1053592-55960.fits",
    ],
    ids=["apstar_dr16", "asstar_dr17", "apstar_dr13"],
)
def test_apStar_loader(tmp_path, url):
    with set_temp_cache(path=str(tmp_path)):
        filename = download_file(url, cache=True)
        spectrum = apStar_loader(filename)  # noqa


@pytest.mark.remote_data
@pytest.mark.parametrize(
    "url",
    [
        "https://data.sdss.org/sas/dr16/apogee/spectro/redux/r12/visit/apo25m/N7789/5094/55874/apVisit-r12-5094-55874-123.fits",
        "https://data.sdss.org/sas/dr17/apogee/spectro/redux/dr17/visit/lco25m/293-37-C/11992/59187/asVisit-dr17-11992-59187-227.fits",
        "https://data.sdss.org/sas/dr13/apogee/spectro/redux/r6/apo25m/4925/55702/apVisit-r6-4925-55702-003.fits",
    ],
    ids=["apvisit_dr16", "asvisit_dr17", "apvisit_dr13"],
)
def test_apVisit_loader(tmp_path, url):
    with set_temp_cache(path=str(tmp_path)):
        filename = download_file(url, cache=True)
        spectrum = apVisit_loader(filename)  # noqa


@pytest.mark.remote_data
@pytest.mark.parametrize(
    "url",
    [
        "https://data.sdss.org/sas/dr16/apogee/spectro/aspcap/r12/l33/apo25m/N7789/aspcapStar-r12-2M00005414+5522241.fits",
        "https://data.sdss.org/sas/dr17/apogee/spectro/aspcap/dr17/synspec_rev1/lco25m/293-37-C/aspcapStar-dr17-2M03174068-7748459.fits",
    ],
    ids=["aspcap_dr16_apo", "aspcap_dr17_lco"],
)
def test_aspcapStar_loader(tmp_path, url):
    with set_temp_cache(path=str(tmp_path)):
        filename = download_file(url, cache=True)
        spectrum = aspcapStar_loader(filename)  # noqa
