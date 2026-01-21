from functools import lru_cache
from typing import Callable, Optional

from astropy.nddata import NDIOMixin

__all__ = ['SpectrumList']


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
        self._id_map: Optional[dict[str, int]] = None

        self._lazy_loader: Optional[Callable] = None
        self._lazy_cache_size: int = 0
        self._lazy_labels: Optional[list[str]] = None

    @property
    def is_lazy(self) -> bool:
        return self._lazy_loader is not None

    @property
    def n_loaded(self) -> int:
        if not self.is_lazy:
            return len(self)

        return sum(1 for x in super().__iter__() if x is not _placeholder)

    def set_id_map(self, id_map: dict[str, int]):
        self._id_map = dict(id_map)

    def _resolve_key(self, key: str) -> int:
        if key.isdigit() and (self._id_map is None or key not in self._id_map):
            return int(key)

        # Otherwise, it must be provided by the mapping.
        if self._id_map is None:
            raise KeyError("No id mapping provided for alternate indexing")

        if key not in self._id_map:
            raise KeyError(f"Key '{key}' not found in id mapping, and cannot resolve to a list index.")

        return self._id_map[key]

    def __getitem__(self, value):
        if isinstance(value, str):
            value = self._resolve_key(value)

        # Preserve normal list slice behavior (and avoid int(slice) errors)
        if isinstance(value, slice):
            return SpectrumList([self[i] for i in range(*value.indices(len(self)))])

        if self.is_lazy:
            return self._lazy_get(int(value))

        return super().__getitem__(value)

    def _lazy_get(self, idx: int):
        if idx < 0:
            idx = len(self) + idx
        if idx < 0 or idx >= len(self):
            raise IndexError("list index out of range")

        current = super().__getitem__(idx)
        if current is not _placeholder:
            return current

        val = self._lazy_loader(idx)
        self[idx] = val

        return val

    def _lazy_repr(self, ii: int):
        """Return repr for item i without triggering lazy loading."""
        item = super().__getitem__(ii)

        if item is not _placeholder:
            return repr(item)

        labels = self._lazy_labels
        if labels is not None and ii < len(labels):
            return repr(labels[ii])

        return repr(item)

    def __repr__(self) -> str:
        if not self.is_lazy:
            return super().__repr__()

        return f"[{', '.join(self._lazy_repr(i) for i in range(len(self)))}]"

    @classmethod
    def from_lazy(cls, *, length: int, loader: Callable, cache_size: int = 0, labels: list = None) -> "SpectrumList":
        """
        Create a lazy-loading `SpectrumList`.

        Parameters
        ----------
        length : int
            The number of spectra in the list.
        loader : Callable
            A function that takes a single integer index and returns the
            corresponding `~specutils.Spectrum` object.
        cache_size : int, optional
            The number of spectra to cache in memory, by default 10.

        Returns
        -------
        SpectrumList
            A lazy-loading `SpectrumList`.
        """
        speclist = cls([_placeholder] * length)
        cache_size = max(int(cache_size), 0)
        if cache_size > 0:
            loader = lru_cache(maxsize=cache_size)(loader)
        speclist._lazy_loader = loader

        speclist._lazy_labels = labels
        return speclist
