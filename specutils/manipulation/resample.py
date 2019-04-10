from abc import ABC, abstractmethod

import numpy as np
from astropy.units import Quantity

from ..spectra import Spectrum1D

class ResampleBase(ABC):
    """
    Base class for resample classes.  The algorithms and needs for difference
    resamples will vary quite a bit, so this class is relatively sparse.
    """

    @abstractmethod
    def __call__(self, orig_spectrum, fin_lamb):
        """
        Return the resulting `~specutils.Spectrum1D` of the resampling.
        """
        return NotImplemented

    @abstractmethod
    def resample(self, orig_spectrum, fin_lamb, weights='unc'):
        """
        Workhorse method that will return the resampled Spectrum1D
        object.
        """
        return NotImplemented

    @staticmethod
    def bin_edges(x):
        """
        Calculate the bin edge values of an input dispersion axis. Input values
        are assumed to be the center of the bins.

        Parameters
        ----------
        x : ndarray
            The input dispersion axis values.

        Returns
        -------
        edges : ndarray
            Calcualated bin edges, including left and right most bin edges.
        """
        inside_edges = (x[1:] + x[:-1]) / 2
        edges = np.insert(inside_edges, 0, 2 * x[0] - inside_edges[0])
        edges = np.append(edges, 2 * x[-1] - inside_edges[-1])

        return edges


class FluxConservingResample(ResampleBase):
    """
    todo: fill out this docstring
    """

    def __call__(self, orig_spectrum, fin_lamb, weights='unc'):
        """
        Return the resulting `~specutils.Spectrum1D` of the resampling.
        """
        return self.resample(orig_spectrum, fin_lamb, weights=weights)

    def _resample_matrix(self, orig_lamb, fin_lamb, weights=None):
        """
        Create a re-sampling matrix to be used in re-sampling spectra in a way
        that conserves flux. This is adapted from *this* paper (ref TBD). This
        code was heavily influenced by Nick Earl's resample rough draft

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
        # Lower bin and upper bin edges
        orig_edges = self.bin_edges(orig_lamb)
        fin_edges = self.bin_edges(fin_lamb)

        # I could get rid of these alias variables,
        # but it does add readability
        orig_low = orig_edges[:-1]
        fin_low = fin_edges[:-1]
        orig_upp = orig_edges[1:]
        fin_upp = fin_edges[1:]

        # Here's the real work in figuring out the bin overlaps
        # i.e., contribution of each original bin to the resampled bin
        l_inf = np.where(orig_low > fin_low[:, np.newaxis],
                         orig_low, fin_low[:, np.newaxis])
        l_sup = np.where(orig_upp < fin_upp[:, np.newaxis],
                         orig_upp, fin_upp[:, np.newaxis])

        # If no weights specified, use all ones
        if not weights:
            weights = np.ones(len(orig_lamb))

        resamp_mat = (l_sup - l_inf).clip(0)
        resamp_mat *= (orig_upp - orig_low)

        # set bins that don't overlap 100% with original bins
        # to zero by checking edges, and applying generated mask
        left_clip = np.where(fin_edges[:-1] - orig_edges[0] < 0, 0, 1)
        right_clip = np.where(orig_edges[-1] - fin_edges[1:] < 0, 0, 1)
        keep_overlapping_matrix = left_clip * right_clip

        resamp_mat *= keep_overlapping_matrix[:, np.newaxis]
        resamp_mat *= weights

        return resamp_mat

    def resample(self, orig_spectrum, fin_lamb,  weights='unc'):
        """
        Create a re-sampling matrix to be used in re-sampling spectra in a way
        that conserves flux. This is adapted from *this* paper (ref TBD)

        Parameters
        ----------
        orig_spectrum : `~specutils.Spectrum1D`
            The original 1D spectrum.
        fin_lamb : ndarray
            The desired dispersion array.
        weights : str or ndarray
            Uncertainty handeling.  The default string 'unc' will use the
            inverse variance from the uncertainty as the weights. No uncertainty
            will be calculated if set to None.  If ndarray is provided, this
            will be used for weights.

        Returns
        -------
        resample_spectrum : `~specutils.Spectrum1D`
            An output spectrum containing the resampled `~specutils.Spectrum1D`
        """

        # Check if units on original spectrum and new wavelength (if defined)
        # match
        if isinstance(fin_lamb, Quantity):
            if orig_spectrum.wavelength.unit != fin_lamb.unit:
                return ValueError("Original spectrum dispersion grid and new"
                                  "dispersion grid must have the same units.")

        # todo: Generate uncertainty if none provided and weights='unc' chosen?
        # todo: We should populate an output uncertainty, but what do we do if input
        # spectrum has populated uncertainty, but a different weighting is
        # provided in the weights keyword?

        # If weight is provided, turn into float array and check if it's the
        # same length as input spectrum
        if weights not in (None,'unc'):
            try:
                weights = np.asfarray(weights)
            except ValueError as e:
                return ValueError("Could not convert weight array to Ndarray: {}".format(e))

            if len(weights) != len(orig_spectrum.flux):
                return ValueError("Provided weight array must be the same "
                                  "length as the input spectrum.")

        # todo, line 74 doesn't like the inputs being quantity objects, may
        # want to look into this more in the future
        resample_grid = self._resample_matrix(np.array(orig_spectrum.wavelength),
                                              np.array(fin_lamb), weights=weights)
        out_flux = np.sum(orig_spectrum.flux * resample_grid, axis=1) / np.sum(
            resample_grid, axis=1)

        # todo: for now, use the units from the pre-resampled
        # spectra, although if a unit is defined for fin_lamb and it doesn't
        # match the input spectrum it won't work right, will have to think
        # more about how to handle that... could convert before and after
        # calculation, which is probably easiest. Matrix math algorithm is
        # geometry based, so won't work to just let quantity math handle it.
        resampled_spectrum = Spectrum1D(np.nan_to_num(out_flux),
                                        np.array(fin_lamb) * orig_spectrum.wavelength.unit)

        return resampled_spectrum