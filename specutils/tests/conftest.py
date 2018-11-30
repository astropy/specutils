import os
import time
import urllib

import pytest
import tempfile


NUM_ATTEMPTS = 3

remote_access = lambda argvals: pytest.mark.parametrize(
    'remote_data_path', argvals, indirect=True, scope='function')


@pytest.fixture(scope='module')
def remote_data_path(request):
    """
    Remotely access the Zenodo deposition archive to retrieve the versioned
    test data.
    """
    # Make use of configuration option from pytest-remotedata in order to
    # control access to remote data.
    if request.config.getoption('remote_data', 'any') != 'any':
        pytest.skip()

    file_id, file_name = request.param.values()

    url = "https://zenodo.org/record/{}/files/{}?download=1".format(
        file_id, file_name)

    for _ in range(NUM_ATTEMPTS):
        with urllib.request.urlopen(url) as r:
            # Take a brief timeout if the request failed
            if r.status != 200:
                time.sleep(3)
                continue

            # Create a temporary directory that is automatically cleaned up
            # when the context is exited, removing any temporarily stored data.
            with tempfile.TemporaryDirectory() as tmp_dir:
                file_path = os.path.join(tmp_dir, file_name)
                with open(file_path, 'wb') as tmp_file:
                    tmp_file.write(r.read())
                    yield file_path
                    break

    else:
        raise ConnectionError("Failed to access remote data at {}".format(url))
