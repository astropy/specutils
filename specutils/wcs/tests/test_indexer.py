from ..specwcs import Indexer
import numpy as np


def test_length():
    ind = Indexer(slice(2, 5, 2))
    assert ind.length == 2
    ind = Indexer(slice(0, 4, 1))
    assert ind.length == 4
    ind = Indexer(slice(0, 4))
    assert ind.length == 4
    ind = Indexer(slice(5, 0, -1))
    assert ind.length == 5
    ind = Indexer(slice(2, 6, 2))
    assert ind.length == 2
    ind = Indexer(slice(6, 2, -2))
    assert ind.length == 2
    ind = Indexer(slice(10, 10, 1))
    assert ind.length == 0
    ind = Indexer(slice(4, 25, -1))
    assert ind.length == 0
    ind = Indexer()
    assert ind.length is None


def test_init_slice():
    x = [1, 2, 3, 4, 5, 6, 7]
    ind = Indexer(slice(2, 5, 2))
    np.testing.assert_allclose(x[2:5:2], x.__getitem__(ind.slice))
    ind = Indexer(slice(None, None, -1))
    np.testing.assert_allclose(x[::-1], x.__getitem__(ind.slice))
    ind = Indexer(slice(5, None, -1))
    np.testing.assert_allclose(x[5::-1], x.__getitem__(ind.slice))
    ind = Indexer(slice(10, None, -1))
    np.testing.assert_allclose(x[10::-1], x.__getitem__(ind.slice))
    ind = Indexer(slice(5, 2, -2))
    np.testing.assert_allclose(x[5:2:-2], x.__getitem__(ind.slice))


def test_apply_slice():
    x = [1, 2, 3, 4, 5, 6, 7]
    ind = Indexer()
    np.testing.assert_allclose(x, x.__getitem__(ind.slice))
    ind = Indexer(slice(2, 5, 2))
    ind.apply_slice(slice(0, 2))
    np.testing.assert_allclose(x[2:5:2], x.__getitem__(ind.slice))

