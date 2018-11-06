import pytest
import requests


REMOTE_DATA_NAMES = [

]


@pytest.fixture(scope='module', params=REMOTE_DATA_NAMES)
def remote_data(request):
