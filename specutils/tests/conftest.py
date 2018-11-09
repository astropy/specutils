import pytest
import requests
import os

import tempfile


@pytest.fixture(scope='module')
def remote_data_path(request):
    """
    Remotely access the Zenodo deposition archive to retrieve the versioned
    test data.
    """
    file_id, file_name = request.param.values()

    f = requests.get("https://zenodo.org/record/{}/files/{}?download=1".format(
        file_id, file_name))

    # Create a temporary directory that is automatically cleaned up when the
    # context is exited, removing any temporarily stored data.
    with tempfile.TemporaryDirectory() as tmp_dir:
        file_path = os.path.join(tmp_dir, file_name)

        with open(file_path, 'wb') as tmp_file:
            tmp_file.write(f.content)

            yield file_path