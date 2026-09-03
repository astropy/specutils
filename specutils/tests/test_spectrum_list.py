import string
from specutils.spectra import SpectrumList


labels = list(string.ascii_lowercase[:10])


def test_nonlazy_spectrum_list():
    """test non-lazy SpectrumList behave like a normal list."""
    sl = SpectrumList(range(10))

    assert not sl.is_lazy
    assert len(sl) == 10
    assert sl.n_loaded == 10
    assert sl.labels is None
    assert '1, 2, 3' in repr(sl)


def test_nonlazy_labels():
    """test non-lazy lists can have labels."""
    sl = SpectrumList(range(10))
    labels = {f"item{i}": i for i in range(10)}
    sl.set_id_map(labels)

    assert not sl.is_lazy
    assert len(sl) == 10
    assert sl.n_loaded == 10
    assert sl.labels == labels
    assert sl.labels["item3"] == 3
    # it does not change the repr
    assert 'item1, item2, item3' not in repr(sl)


def test_lazy_spectrum_list_no_labels():
    """test lazy loading items with no labels."""
    # list of items
    items = list(range(10))

    # define an item loader
    def loader(i):
        return items[i]

    sl = SpectrumList.from_lazy(length=len(items), loader=loader)

    # check if lazy
    assert sl.is_lazy
    assert len(sl) == 10
    assert sl.n_loaded == 0
    assert sl.labels is None

    # load an item
    assert sl[3] == 3
    assert sl.n_loaded == 1
    assert 'lazy list: 1 items loaded' in repr(sl)

    # check the first item isn't loaded yet
    assert "load a spectrum:\n[<object object at " in repr(sl)

    # load another item
    assert sl[-1] == 9
    assert sl.n_loaded == 2
    assert 'lazy list: 2 items loaded' in repr(sl)


def test_lazy_spectrum_list_with_labels():
    """test lazy loading items with labels."""

    # set list of items and labels
    items = list(range(10))
    labels = [f"item{i}" for i in items]

    # define an item loader
    def loader(i):
        return items[i]

    sl = SpectrumList.from_lazy(length=len(items), loader=loader, labels=labels)

    assert sl.is_lazy
    assert sl.n_loaded == 0
    assert sl.labels == labels
    assert isinstance(sl.labels, list)
    assert "'item1', 'item2', 'item3'" in repr(sl)

    assert sl["item3"] == 3
    assert sl.n_loaded == 1
    assert "'item1', 'item2', 3" in repr(sl)
