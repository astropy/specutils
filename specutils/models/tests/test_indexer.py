from specutils.models.Indexer import Indexer
import numpy as np
import astropy.io.ascii as ascii
from astropy import units as u
import os
from specutils.io import read_fits


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_spectrum_slicing():
    iraf = ascii.read(data_path('gbt_1d_iraf_read.dat'), names=['wave', 'flux'])
    spec = read_fits.read_fits_spectrum1d(data_path('gbt_1d.fits'))
    np.testing.assert_allclose(spec.dispersion.value, iraf['wave'])
    np.testing.assert_allclose(spec.flux, iraf['flux'])
    sliced = spec.slice_index(stop=100)
    np.testing.assert_allclose(sliced.dispersion.value, iraf["wave"][:100])
    np.testing.assert_allclose(sliced.flux, iraf["flux"][:100])
    np.testing.assert_allclose(spec.dispersion.value, iraf['wave'])
    np.testing.assert_allclose(spec.flux, iraf['flux'])
    sliced = sliced.slice_index(start=-1, step=-1)
    np.testing.assert_allclose(sliced.dispersion.value, iraf["wave"][99::-1])
    np.testing.assert_allclose(sliced.flux, iraf["flux"][99::-1])


def test_length():
    ind = Indexer(2, 5, 2)
    assert ind.length == 2
    ind = Indexer(0, 4, 1)
    assert ind.length == 4
    ind = Indexer(0, 4)
    assert ind.length == 4
    ind = Indexer(5, -1, -1)
    assert ind.length == 6
    ind = Indexer(2, 6, 2)
    assert ind.length == 2
    ind = Indexer(6, 2, -2)
    assert ind.length == 2
    ind = Indexer(10, 10, 1)
    assert ind.length == 0
    ind = Indexer(4, 25, -1)
    assert ind.length == 0


def test_init_slice():
    x = np.arange(101)
    ind = Indexer(2, 5, 2)
    np.testing.assert_allclose(x[2:5:2], x[ind()])
    ind = Indexer(0, 100)
    np.testing.assert_allclose(x[0:100], x[ind()])
    ind = Indexer(100, -1, -1)
    np.testing.assert_allclose(x[100::-1], x[ind()])
    ind = Indexer(10, 3, -2)
    np.testing.assert_allclose(x[10:3:-2], x[ind()])
    ind = Indexer(5, 2, -2)
    np.testing.assert_allclose(x[5:2:-2], x[ind()])


def test_apply_slice():
    x = np.arange(1024)
    ind = Indexer(2, 5, 2)
    ind2 = ind[0:2]
    np.testing.assert_allclose(x[2:5:2], x[ind()])
    np.testing.assert_allclose(x[2:5:2], x[ind2()])
    ind = Indexer(0, 1024)
    ind2 = ind[::-1]
    np.testing.assert_allclose(x, x[ind()])
    np.testing.assert_allclose(x[::-1], x[ind2()])
    np.testing.assert_allclose(x[100::-1], x[ind[100::-1]()])
    np.testing.assert_allclose(x[:100:-1], x[ind[:100:-1]()])
    np.testing.assert_allclose(x[200:100:-1], x[ind[200:100:-1]()])
    np.testing.assert_allclose(x[242:100:-5], x[ind[242:100:-5]()])
    np.testing.assert_allclose(x[242:100], x[ind[242:100]()])
    np.testing.assert_allclose(x[:100], x[ind[:100]()])
    np.testing.assert_allclose(x[:100:3], x[ind[:100:3]()])
    np.testing.assert_allclose(x[:], x[ind[:]()])
    np.testing.assert_allclose(x[50:], x[ind[50:]()])
    np.testing.assert_allclose(x[50::2], x[ind[50::2]()])
    np.testing.assert_allclose(x[:-1:-1], x[ind[:-1:-1]()])
    np.testing.assert_allclose(x[:-100:-1], x[ind[:-100:-1]()])
    np.testing.assert_allclose(x[-200:-100], x[ind[-200:-100]()])
    np.testing.assert_allclose(x[-200:-100:-1], x[ind[-200:-100:-1]()])
    np.testing.assert_allclose(x[-100:-200:-1], x[ind[-100:-200:-1]()])
    np.testing.assert_allclose(x[-100::-1], x[ind[-100::-1]()])
    np.testing.assert_allclose(x[-100:], x[ind[-100:]()])
    np.testing.assert_allclose(x[-100::], x[ind[-100::]()])

