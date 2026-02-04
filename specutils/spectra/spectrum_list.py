from functools import lru_cache
from typing import Callable, Optional

from astropy.nddata import NDIOMixin

__all__ = ['SpectrumList']

# a temporary placeholder object for lists
_placeholder = object()


class SpectrumList(list, NDIOMixin):
    """
    A list that is used to hold a list of `~specutils.Spectrum` objects

    The primary purpose of this class is to allow loaders to return a list of
    spectra that have different shapes.  For spectra that have the same shape
    but different spectral axes, see `~specutils.SpectrumCollection`.  For
    a spectrum or spectra that all share the same spectral axis, use
    `~specutils.Spectrum`.  For more on this topic, see
    :ref:`specutils-representation-overview`.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Mapping of alternate string ids to list index
        self._id_map: Optional[dict[str, int]] = None

        # Parameters for lazy loading
        self._lazy_loader: Optional[Callable] = None
        self._lazy_cache_size: int = 0
        self._lazy_labels: Optional[list[str]] = None

    @property
    def is_lazy(self) -> bool:
        """Whether the SpectrumList is in lazy-loading mode"""
        return self._lazy_loader is not None

    @property
    def n_loaded(self) -> int:
        """Number of spectra in the list currently loaded"""
        if not self.is_lazy:
            return len(self)

        return sum(1 for x in super().__iter__() if x is not _placeholder)

    def set_id_map(self, id_map: dict[str, int]):
        """Set a mapping of alternate string labels to list indices.

        This allows accessing items in the list using string labels, e.g.
        source ids, in addition to int indices.

        Parameters
        ----------
        id_map : dict[str, int]
            Mapping of string keys to list indices
        """
        self._id_map = dict(id_map)

    def _resolve_key(self, key: str) -> int:
        """Resolve a string key to a list index"""

        # return normal list index
        if key.isdigit() and (self._id_map is None or key not in self._id_map):
            return int(key)

        # otherwise it must be provided by the mapping
        if self._id_map is None:
            raise KeyError("No id mapping provided for alternate indexing")

        if key not in self._id_map:
            raise KeyError(f"Key '{key}' not found in id mapping, and cannot resolve to a list index.")

        return self._id_map[key]

    def __getitem__(self, value: str | int):
        """Retrieve items from the list normally or lazily"""

        # resolve any string key
        if isinstance(value, str):
            value = self._resolve_key(value)

        # preserve original slice behaviour
        if isinstance(value, slice):
            return SpectrumList([self[i] for i in range(*value.indices(len(self)))])

        # use lazy item getter
        if self.is_lazy:
            return self._lazy_get(int(value))

        # use normal item getter
        return super().__getitem__(value)

    def _lazy_get(self, idx: int):
        """Lazily retrieve an item from the list"""

        # preserve original reverse slicing
        if idx < 0:
            idx = len(self) + idx
        if idx < 0 or idx >= len(self):
            raise IndexError("list index out of range")

        # get the current item and check if it's a placeholder object
        # if not, then return it
        current = super().__getitem__(idx)
        if current is not _placeholder:
            return current

        # use the lazy loader to get the spectrum object
        # replace the placeholder item with it
        val = self._lazy_loader(idx)
        self[idx] = val

        return val

    def __repr__(self) -> str:
        """Build string repr the list

        Non-lazy lists have normal reprs.  Lazy lists use placeholder values,
        or optional labels if provided.  Once the item is loaded, the normal
        item repr is used.
        """
        # use normal repr
        if not self.is_lazy or self.n_loaded == len(self):
            return super().__repr__()

        # build the lazy repr
        prefix = f"lazy list: {self.n_loaded} items loaded; access an index to load a spectrum:\n"
        return prefix + f"[{', '.join(self._lazy_repr(i) for i in range(len(self)))}]"

    def _lazy_repr(self, ii: int):
        """Return item repr without triggering a lazy load"""

        # get item
        item = super().__getitem__(ii)

        # use normal item repr
        if item is not _placeholder:
            return repr(item)

        # use placeholder or optional label repr
        labels = self._lazy_labels
        if labels is not None and ii < len(labels):
            return repr(labels[ii])

        return repr(item)

    @classmethod
    def from_lazy(
        cls, length: int, loader: Callable, cache_size: Optional[int] = None, labels: list = None
    ) -> "SpectrumList":
        """Construct a lazy-loading SpectrumList.

        Constructs a Spectrumlist using placeholder objects and sets a
        loader callable used to instantiate Spectrum objects lazily on item
        get. Once a Spectrum is loaded, it replaces the placeholder item, and
        repeated access does not re-run the loader.

        If cache_size is specified, then the loader is also cached with
        an lru_cache.

        Parameters
        ----------
        length : int
            Total number of spectra in the list.
        loader : Callable
            Callable taking an int index and returns a Spectrum.
        cache_size : int or None, optional
            If provided, wraps the loader in an lru_cache of this
            size.
        labels : list, optional
            Optional list of placeholder display labels shown by the repr
            before spectra are materialized.

        Returns
        -------
        SpectrumList
            A lazy-loadable SpectrumList

        """
        # create placeholder list
        speclist = cls([_placeholder] * length)

        # optionally cache the loader
        if cache_size:
            cache_size = max(int(cache_size), 0)
            loader = lru_cache(maxsize=cache_size)(loader)

        # set lazy parameters
        speclist._lazy_loader = loader
        speclist._lazy_labels = labels
        return speclist
