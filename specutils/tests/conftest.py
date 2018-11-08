import pytest
import requests
import os

import tempfile


REMOTE_DATA_IDS = {
    'cos_fuv': ('1481119', 'COS_FUV.fits')
}

REMOTE_DATA_NAMES = {
    'COS_FUV.fits': '1481119'
}

ACCESS_TOKEN = os.environ.get("ACCESS_TOKEN")


@pytest.fixture(scope='module', params=REMOTE_DATA_IDS.values(), ids=REMOTE_DATA_IDS.keys())
def remote_data(request):
    if ACCESS_TOKEN is None:
        raise Exception("Not access token set in the environment.")

    r = requests.get("https://zenodo.org/api/deposit/depositions/{}/files/{}".format(
        request.param[0], request.param[1]), params={'access_token': ACCESS_TOKEN})

    f = requests.get(r.json()['links']['download'], params={'access_token': ACCESS_TOKEN})

    with tempfile.TemporaryDirectory() as tmp_dir:
        file_path = os.path.join(tmp_dir, request.param[1])
        open(file_path, 'wb').write(f.content)

        yield file_path