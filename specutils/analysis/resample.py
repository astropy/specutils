import logging
import os

import numpy as np
import astropy.units as u

__all__ = ['Resample']


class Resample(object):
    """
    Resample model which can be used with compound model objects.
    """

    u.quantity_input(new_dispersion=['speed', 'length'])
    def __init__(self, new_dispersion):
        self._new_dispersion = new_dispersion
        self._disp_unit = self._new_dispersion.unit

    def __call__(self, x, y):
        _remat = resample(x, self._new_dispersion)

        return np.dot(_remat, y)


def resample(orig_lamb, fin_lamb, unit=None, force=None, **kwargs):
    """
    Abstraction function to decide whether to use uniform or non-uniform
    resampling functions.

    Parameters
    ----------
    orig_lamb : ndarray
        Original dispersion grid.
    fin_lamb : ndarray
        Final dispersion grid.
    kwargs : dict
        Extra keyword arguments to pass to functions.

    Returns
    -------
    resample_mat : ndarray
        The resampling matrix to be applied to data arrays.
    """

    orig_lamb = u.Quantity(orig_lamb, unit=unit)
    fin_lamb = u.Quantity(fin_lamb, unit=unit)

    orig_space = (orig_lamb[1:] - orig_lamb[:-1]).value
    fin_space = (fin_lamb[1:] - fin_lamb[:-1]).value

    if np.allclose(orig_space, orig_space[0]) and np.allclose(fin_space, fin_space[0]):
        logging.info("Re-sampling: original and final grids are uniform.")
        mat = _uniform_matrix(orig_lamb, fin_lamb)
    else:
        logging.info("Re-sampling: original and final grids are non-uniform.")
        mat = _nonuniform_matrix(orig_lamb, fin_lamb, **kwargs)

    return mat


def _uniform_matrix(orig_lamb, fin_lamb):
    """
    Create a re-sampling matrix to be used in re-sampling spectra in a way
    that conserves flux. This is adapted from code created by the SEAGal
    Group.
    .. note:: This method assumes uniform grids.
    Parameters
    ----------
    orig_lamb : ndarray
        The original dispersion array.
    fin_lamb : ndarray
        The desired dispersion array.
    Returns
    -------
    resample_mat : ndarray
        An [[N_{fin_lamb}, M_{orig_lamb}]] matrix.
    """

    # Get step size
    delta_orig = orig_lamb[1] - orig_lamb[0]
    delta_fin = fin_lamb[1] - fin_lamb[0]

    n_orig_lamb = len(orig_lamb)
    n_fin_lamb = len(fin_lamb)

    # Lower bin and upper bin edges
    orig_low = orig_lamb - delta_orig * 0.5
    orig_upp = orig_lamb + delta_orig * 0.5
    fin_low = fin_lamb - delta_fin * 0.5
    fin_upp = fin_lamb + delta_fin * 0.5

    # Create re-sampling matrix
    resamp_mat = np.zeros(shape=(n_fin_lamb, n_orig_lamb))

    for i in range(n_fin_lamb):
        # Calculate the contribution of each original bin to the
        # resampled bin
        l_inf = np.where(orig_low > fin_low[i], orig_low, fin_low[i])
        l_sup = np.where(orig_upp < fin_upp[i], orig_upp, fin_upp[i])

        # Interval overlap of each original bin for current resampled
        # bin; negatives clipped
        dl = (l_sup - l_inf).clip(0)

        # This will only happen at the edges of lorig.
        # Discard resampled bin if it's not fully covered (> 99%) by the
        #  original bin -- only happens at the edges of the original bins
        if 0 < dl.sum() < 0.99 * delta_fin.value:
            dl = 0 * orig_lamb

        resamp_mat[i, :] = dl

    resamp_mat = resamp_mat / delta_fin

    return resamp_mat


def _nonuniform_matrix(orig_lamb, fin_lamb, extrapolate=False):
    """
    Compute re-sampling matrix R_o2r, useful to convert a spectrum sampled at
    wavelengths lorig to a new grid lresamp. Here, there is no necessity to
    have constant gris as on :func:`_uniform_matrix`. This is adapted from code
    created by the SEAGal Group.
    .. warning:: lorig and lresam MUST be on ascending order!
    Parameters
    ----------
    orig_lamb : array_like
        Original spectrum lambda array.
    fin_lamb : array_like
        Spectrum lambda array in which the spectrum should be sampled.
    extrapolate : boolean, optional
        Extrapolate values, i.e. values for fin_lamb < orig_lamb[0] are set to
        match orig_lamb[0] and values for fin_lamb > orig_lamb[-1] are set to
        match orig_lamb[-1].
    Returns
    -------
    resample_mat : ndarray
        Resample matrix.
    """
    matrix = np.zeros(shape=(len(fin_lamb), len(orig_lamb)))

    # Define lambda ranges (low, upp) for original and resampled.
    lo_low = np.zeros(len(orig_lamb))
    lo_low[1:] = (orig_lamb[1:] + orig_lamb[:-1]) / 2
    lo_low[0] = orig_lamb[0] - (orig_lamb[1] - orig_lamb[0]) / 2

    lo_upp = np.zeros(len(orig_lamb))
    lo_upp[:-1] = lo_low[1:]
    lo_upp[-1] = orig_lamb[-1] + (orig_lamb[-1] - orig_lamb[-2]) / 2

    lr_low = np.zeros(len(fin_lamb))
    lr_low[1:] = (fin_lamb[1:] + fin_lamb[:-1]) / 2
    lr_low[0] = fin_lamb[0] - (fin_lamb[1] - fin_lamb[0]) / 2

    lr_upp = np.zeros(len(fin_lamb))
    lr_upp[:-1] = lr_low[1:]
    lr_upp[-1] = fin_lamb[-1] + (fin_lamb[-1] - fin_lamb[-2]) / 2

    # Iterate over resampled fin_lamb vector
    for i in range(len(fin_lamb)):

        # Find in which bins fin_lamb bin within orig_lamb bin
        bins_resamp = np.where((lr_low[i] < lo_upp) & (lr_upp[i] > lo_low))[0]

        # On these bins, eval fraction of resampled bin is within original bin.
        for j in bins_resamp:
            aux = 0

            d_lr = lr_upp[i] - lr_low[i]
            d_lo = lo_upp[j] - lo_low[j]
            d_ir = lo_upp[j] - lr_low[i]  # common section on the right
            d_il = lr_upp[i] - lo_low[j]  # common section on the left

            # Case 1: resampling window is smaller than or equal to the
            # original window. This is where the bug was: if an original bin
            # is all inside the resampled bin, then all flux should go into it,
            # not then d_lr/d_lo fraction.
            if (lr_low[i] > lo_low[j]) & (lr_upp[i] < lo_upp[j]):
                aux += 1.

            # Case 2: resampling window is larger than the original window.
            if (lr_low[i] < lo_low[j]) & (lr_upp[i] > lo_upp[j]):
                aux += d_lo / d_lr

            # Case 3: resampling window is on the right of the original window.
            if (lr_low[i] > lo_low[j]) & (lr_upp[i] > lo_upp[j]):
                aux += d_ir / d_lr

            # Case 4: resampling window is on the left of the original window.
            if (lr_low[i] < lo_low[j]) & (lr_upp[i] < lo_upp[j]):
                aux += d_il / d_lr

            matrix[i, j] += aux

    # Fix extremes: extrapolate if needed
    if extrapolate:
        bins_extrapl = np.where((lr_low < lo_low[0]))[0]
        bins_extrapr = np.where((lr_upp > lo_upp[-1]))[0]

        if len(bins_extrapl) > 0 and len(bins_extrapr) > 0:
            io_extrapl = np.where((lo_low >= lr_low[bins_extrapl[0]]))[0][0]
            io_extrapr = np.where((lo_upp <= lr_upp[bins_extrapr[0]]))[0][-1]

            matrix[bins_extrapl, io_extrapl] = 1.
            matrix[bins_extrapr, io_extrapr] = 1.

    return matrix
