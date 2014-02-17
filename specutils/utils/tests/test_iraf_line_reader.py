from .. import IRAF_id_reader
import os

def data_path(filename):
        data_dir = os.path.join(os.path.dirname(__file__), 'data')
        return os.path.join(data_dir, filename)

def test_iraf_line_reader():
    IRAF_id_reader.IRAF_identify_reader(data_path('idexample'))
