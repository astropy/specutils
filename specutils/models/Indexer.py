from astropy.modeling import Model
import math
import numpy as np
from astropy.modeling.parameters import Parameter

class Indexer(Model):
    """
    Indexer class provides an alternate slicing mechanism for WCS'es. As WCS'es
    are defined for negative integers as well, this indexer class can be
    initialized with negative indices. Unlike python slice objects, a
    negative index is the same as a positive index. However, this is only true
    for initialization. Once initialized, the indexer can be used like normal
    python slicing, with negative index pointing to length - index.

    Parameters
    --------------
    start: int
        the start point of the indexer, can be any integer
    stop: int
        the point up to which (not including) the index is valid, can be any
        integer
    step: int, optional
        the jump between two indices

    Raises
    -----------
    ValueError
        if step is given as zero
    """
    # start = Parameter()
    # stop = Parameter()
    # step = Parameter()
    # length = Parameter()

    def __init__(self, start, stop, step=1):
        self._n_models = 1
        # self.param_names = ["start", "stop", "step", "length"]
        if step == 0:
            raise ValueError("slice step cannot be zero")
        self.start = start
        self.stop = stop
        self.step = step

    def _parse_slice(self, slice):
        """
        Internal method to extract the start, stop and stop from the given slice
        object
        """
        if slice.step == 0:
            raise ValueError("slice step cannot be zero")
        step = slice.step if slice.step is not None else 1
        if step < 0:
            default_start = self.length - 1
            default_stop = -1
        else:
            default_start = 0
            default_stop = self.length
        start = slice.start if slice.start is not None else default_start
        stop = slice.stop if slice.stop is not None else default_stop
        if slice.start is not None:
            start += self.length if start < 0 else 0
            start = 0 if start < 0 else start
            start = self.length - 1 if start >= self.length else start
        if slice.stop is not None:
            stop += self.length if stop < 0 else 0
            stop = -1 if stop < 0 else stop
            stop = self.length if stop > self.length else stop
        return start, stop, step

    def __getitem__(self, d_slice):
        """
        Applies a slice on top of the current slice object, returning a new
        indexer with the modifications.
        NOTE: This method does not modify the current indexer

        Parameters
        -----------
        d_slice: slice object
            The new slice object to be applied
        """
        d_start, d_stop, d_step = self._parse_slice(d_slice)
        n_stop = self.start + d_stop * self.step
        n_start = self.start + d_start * self.step
        n_step = self.step * d_step
        return Indexer(n_start, n_stop, n_step)

    def __call__(self, indices=None):
        """
        Transforms the input indices to the correct indices represented by the
        indexer, removes those indexes which are not in range (i.e. beyond, and
        including self.stop)

        Parameters
        -----------
        indices: numpy array
            The indices to be shifted
        """
        if indices is None:
            indices = np.arange(self.length)
        return self.__class__.evaluate(indices, self.start,
                                       self.stop, self.step)

    @classmethod
    def evaluate(cls, indices, start, stop, step):
        n_indices = indices * step + start
        if step > 0:
            to_delete = np.where(n_indices >= stop)
        else:
            to_delete = np.where(n_indices <= stop)
        n_indices = np.delete(n_indices, to_delete)
        return n_indices

    @property
    def length(self):
        if (self.start - self.stop) * self.step > 0:
            # the index represents a null list
            return 0
        return int(math.ceil(
            math.fabs((self.start - self.stop) * 1.0 / self.step)))


