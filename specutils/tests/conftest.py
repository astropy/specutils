import pytest
import requests
import os

import tempfile

ZENODO_ACCESS_TOKEN = os.environ.get("ZENODO_ACCESS_TOKEN")


@pytest.fixture(scope='module')
def remote_data_path(request):
    """
    Remotely access the Zenodo deposition archive to retrieve the versioned
    test data.
    """
    if ZENODO_ACCESS_TOKEN is None:
        raise Exception("No access token set in the environment.")

    file_id, file_name = request.param.values()

    # Retrieve information on the deposition files for the given id
    r = requests.get("https://zenodo.org/api/deposit/depositions/{}/files/{}".format(
        file_id, file_name), params={'access_token': ZENODO_ACCESS_TOKEN})

    # Directly download the data file link returned from the deposition info
    f = requests.get(r.json()['links']['download'], params={'access_token': ZENODO_ACCESS_TOKEN})

    # Create a temporary directory that is automatically cleaned up when the
    # context is exited, removing any temporarily stored data.
    with tempfile.TemporaryDirectory() as tmp_dir:
        file_path = os.path.join(tmp_dir, file_name)

        with open(file_path, 'wb') as tmp_file:
            tmp_file.write(f.content)

            yield file_path