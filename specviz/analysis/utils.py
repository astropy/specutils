import numpy as np

from ..core.data import Data


def resample(from_data, to_data, copy=False):
    """
    A function to resample a give `Data` or `Layer` object.

    Parameters
    ----------
    data : :class:`Data` or :class:`Layer`
        Object containing the data to be resampled.
    new_lambda : ndarray
        New wavelength grid onto which the `data` object will be resampled.
    copy : bool
        Copy the `data` object.

    Returns
    -------
    new_data : :class:`Data` or :class:`Layer`
        New data object.
    """
    if from_data.dispersion.size > to_data.dispersion.size:
        remat = resample_matrix(from_data.dispersion.value,
                                to_data.dispersion.value)
        flux = np.dot(remat, from_data.data.value)
    else:
        flux = (to_data.dispersion.value, from_data.dispersion.value,
                from_data.data.value)

    new_data = from_data._from_self(flux)

    return new_data


def resample_matrix(self, orig_lamb, fin_lamb):
    """
    Create a resampling matrix to be used in resampling spectra in a way
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
    resample_map : ndarray
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

    # Create resampling matrix
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
        if 0 < dl.sum() < 0.99 * delta_fin:
            dl = 0 * orig_lamb

        resamp_mat[i, :] = dl

    resamp_mat /= delta_fin

    return resamp_mat