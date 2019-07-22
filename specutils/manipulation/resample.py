from abc import ABC, abstractmethod

import numpy as np
from scipy.interpolate import CubicSpline
from astropy.units import Quantity
from astropy.nddata import StdDevUncertainty, VarianceUncertainty, InverseVariance

from ..spectra import Spectrum1D

__all__ = ['ResamplerBase', 'FluxConservingResampler',
           'LinearInterpolatedResampler', 'SplineInterpolatedResampler']


class ResamplerBase(ABC):
    """
    Base class for resample classes.  The algorithms and needs for difference
    resamples will vary quite a bit, so this class is relatively sparse.

    Init paramtere here is not yet hooked up to the rest of the code, but to
    show how we will want to use it in the future.
    """
    def __init__(self, bin_edges='nan_fill'):
        self.bin_edges = bin_edges

    @abstractmethod
    def __call__(self, orig_spectrum, fin_lamb):
        """
        Return the resulting `~specutils.Spectrum1D` of the resampling.
        """
        return NotImplemented

    @abstractmethod
    def resample1d(self, orig_spectrum, fin_lamb):
        """
        Workhorse method that will return the resampled Spectrum1D
        object.
        """
        return NotImplemented

    @staticmethod
    def _calc_bin_edges(x):
        """
        Calculate the bin edge values of an input dispersion axis. Input values
        are assumed to be the center of the bins.

        todo: this should live in the main spectrum object, but we're still
        figuring out the details to that implementation, so leaving here
        for now.

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


class FluxConservingResampler(ResamplerBase):
    """
    This resampling algorithm conserves overall integrated flux (as opposed to
    flux density).
    Algorithm based on the equations documented in the following paper:
    https://ui.adsabs.harvard.edu/abs/2017arXiv170505165C/abstract

    Examples
    --------

    To resample an input spectrum to a user specified dispersion grid using
    a flux conserving algorithm:

    >>> import numpy as np
    >>> import astropy.units as u
    >>> from specutils import Spectrum1D
    >>> from specutils.manipulation import FluxConservingResampler
    >>> input_spectra = Spectrum1D(
    ...     flux=np.array([1, 3, 7, 6, 20]) * u.mJy,
    ...     spectral_axis=np.array([2, 4, 12, 16, 20]) * u.nm)
    >>> resample_grid = np.array([1, 5, 9, 13, 14, 17, 21, 22, 23])
    >>> fluxc_resample = FluxConservingResampler()
    >>> output_spectrum1D = fluxc_resample(input_spectra, resample_grid) # doctest: +IGNORE_OUTPUT

    """

    def __call__(self, orig_spectrum, fin_lamb):
        """
        Return the resulting `~specutils.Spectrum1D` of the resampling.
        """
        return self.resample1d(orig_spectrum, fin_lamb)

    def _resample_matrix(self, orig_lamb, fin_lamb):
        """
        Create a re-sampling matrix to be used in re-sampling spectra in a way
        that conserves flux. This code was heavily influenced by Nick Earl's
        resample rough draft: nmearl@0ff6ef1.

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
        orig_edges = self._calc_bin_edges(orig_lamb)
        fin_edges = self._calc_bin_edges(fin_lamb)

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

        resamp_mat = (l_sup - l_inf).clip(0)
        resamp_mat *= (orig_upp - orig_low)

        # set bins that don't overlap 100% with original bins
        # to zero by checking edges, and applying generated mask
        left_clip = np.where(fin_edges[:-1] - orig_edges[0] < 0, 0, 1)
        right_clip = np.where(orig_edges[-1] - fin_edges[1:] < 0, 0, 1)
        keep_overlapping_matrix = left_clip * right_clip

        resamp_mat *= keep_overlapping_matrix[:, np.newaxis]

        return resamp_mat

    def resample1d(self, orig_spectrum, fin_lamb):
        """
        Create a re-sampling matrix to be used in re-sampling spectra in a way
        that conserves flux. If an uncertainty is present in the input spectra
        it will be propagated through to the final resampled output spectra
        as an InverseVariance uncertainty.

        Parameters
        ----------
        orig_spectrum : `~specutils.Spectrum1D`
            The original 1D spectrum.
        fin_lamb : ndarray
            The desired dispersion array.

        Returns
        -------
        resample_spectrum : `~specutils.Spectrum1D`
            An output spectrum containing the resampled `~specutils.Spectrum1D`
        """

        # Check if units on original spectrum and new wavelength (if defined)
        # match
        if isinstance(fin_lamb, Quantity):
            if orig_spectrum.spectral_axis_unit != fin_lamb.unit:
                return ValueError("Original spectrum dispersion grid and new"
                                  "dispersion grid must have the same units.")

        # todo: Would be good to return uncertainty in type it was provided?
        # todo: add in weighting options

        # Get provided uncertainty into variance
        if orig_spectrum.uncertainty is not None:
            if isinstance(orig_spectrum.uncertainty, StdDevUncertainty):
                pixel_uncer = np.square(orig_spectrum.uncertainty.array)
            elif isinstance(orig_spectrum.uncertainty, VarianceUncertainty):
                pixel_uncer = orig_spectrum.uncertainty.array
            elif isinstance(orig_spectrum.uncertainty, InverseVariance):
                pixel_uncer = np.reciprocal(orig_spectrum.uncertainty.array)
        else:
            pixel_uncer = None

        # todo: Current code doesn't like the inputs being quantity objects, may
        # want to look into this more in the future
        resample_grid = self._resample_matrix(np.array(orig_spectrum.spectral_axis),
                                              np.array(fin_lamb))

        # Now for some broadcasting magic to handle multi dimensional flux inputs
        # Essentially this part is inserting length one dimensions as fillers
        # For example, if we have a (5,6,10) input flux, and an output grid
        # of 3, flux will be broadcast to (5,6,1,10) and resample_grid will
        # Be broadcast to (1,1,3,10).  The sum then reduces down the 10, the
        # original dispersion grid, leaving 3, the new dispersion grid, as
        # the last index.
        new_flux_shape = list(orig_spectrum.flux.shape)
        new_flux_shape.insert(-1, 1)
        in_flux = orig_spectrum.flux.reshape(new_flux_shape)

        ones = [1] * len(orig_spectrum.flux.shape[:-1])
        new_shape_resample_grid = ones + list(resample_grid.shape)
        resample_grid = resample_grid.reshape(new_shape_resample_grid)

        # Calculate final flux
        out_flux = np.sum(in_flux * resample_grid, axis=-1) / np.sum(
            resample_grid, axis=-1)

        # Calculate output uncertainty
        if pixel_uncer is not None:
            pixel_uncer = pixel_uncer.reshape(new_flux_shape)

            out_variance = np.sum(pixel_uncer * resample_grid**2, axis=-1) / np.sum(
                resample_grid**2, axis=-1)
            out_uncertainty = InverseVariance(np.reciprocal(out_variance))
        else:
            out_uncertainty = None

        # todo: for now, use the units from the pre-resampled
        # spectra, although if a unit is defined for fin_lamb and it doesn't
        # match the input spectrum it won't work right, will have to think
        # more about how to handle that... could convert before and after
        # calculation, which is probably easiest. Matrix math algorithm is
        # geometry based, so won't work to just let quantity math handle it.
        resampled_spectrum = Spectrum1D(flux=out_flux,
                                        spectral_axis=np.array(fin_lamb) * orig_spectrum.spectral_axis_unit,
                                        uncertainty=out_uncertainty)

        return resampled_spectrum


class LinearInterpolatedResampler(ResamplerBase):
    """
    Resample a spectrum onto a new ``spectral_axis`` using linear interpolation.

    Examples
    --------

    To resample an input spectrum to a user specified dispersion grid using
    linear interpolation:

    >>> import numpy as np
    >>> import astropy.units as u
    >>> from specutils import Spectrum1D
    >>> from specutils.manipulation import LinearInterpolatedResampler
    >>> input_spectra = Spectrum1D(
    ...     flux=np.array([1, 3, 7, 6, 20]) * u.mJy,
    ...     spectral_axis=np.array([2, 4, 12, 16, 20]) * u.nm)
    >>> resample_grid = np.array([1, 5, 9, 13, 14, 17, 21, 22, 23])
    >>> fluxc_resample = LinearInterpolatedResampler()
    >>> output_spectrum1D = fluxc_resample(input_spectra, resample_grid) # doctest: +IGNORE_OUTPUT
    """
    def __init__(self, bin_edges='nan_fill'):
        super().__init__(bin_edges)

    def __call__(self, orig_spectrum, fin_lamb):
        """
        Return the resulting `~specutils.Spectrum1D` of the resampling.
        """
        return self.resample1d(orig_spectrum, fin_lamb)

    def _interpolation(self, orig_dispersion, flux, fin_lamb):
        """
        Use specified interpolation to calculated resampled
        flux.

        Parameters
        ----------
        orig_dispersion : ndarray
            The original dispersion array.
        flux: ndarray
            The flux array from the input Spectrum1D
        fin_lamb : ndarray
            The desired dispersion array.

        Returns
        -------
        resample_flux : ndarray
            The resampled flux array generated from the interpolation.

        """
        return np.interp(fin_lamb, orig_dispersion, flux, left=np.nan, right=np.nan)

    def resample1d(self, orig_spectrum, fin_lamb):
        """
        Call interpolation, repackage new spectra


        Parameters
        ----------
        orig_spectrum : `~specutils.Spectrum1D`
            The original 1D spectrum.
        fin_lamb : ndarray
            The desired dispersion array.

        Returns
        -------
        resample_spectrum : `~specutils.Spectrum1D`
            An output spectrum containing the resampled `~specutils.Spectrum1D`
        """
        if orig_spectrum.uncertainty is not None:
            warn("Linear interpolation currently does not propogate uncertainties")
        out_flux = self._interpolation(orig_spectrum.spectral_axis, orig_spectrum.flux,
                                  fin_lamb)

        # todo: for now, use the units from the pre-resampled
        # spectra, although if a unit is defined for fin_lamb and it doesn't
        # match the input spectrum it won't work right, will have to think
        # more about how to handle that... could convert before and after
        # calculation, which is probably easiest. Matrix math algorithm is
        # geometry based, so won't work to just let quantity math handle it.
        # todo: handle uncertainties for interpolated cases.
        resampled_spectrum = Spectrum1D(flux=out_flux * orig_spectrum.flux.unit,
                                        spectral_axis=np.array(fin_lamb) * orig_spectrum.spectral_axis_unit)

        return resampled_spectrum


class SplineInterpolatedResampler(ResamplerBase):
    """
    This resample algorithim uses a cubic spline interpolator.  In the future
    this can be expanded to use splines of different degrees.

    Examples
    --------

    To resample an input spectrum to a user specified dispersion grid using
    a cubic spline interpolator:

    >>> import numpy as np
    >>> import astropy.units as u
    >>> from specutils import Spectrum1D
    >>> from specutils.manipulation import SplineInterpolatedResampler
    >>> input_spectra = Spectrum1D(
    ...     flux=np.array([1, 3, 7, 6, 20]) * u.mJy,
    ...     spectral_axis=np.array([2, 4, 12, 16, 20]) * u.nm)
    >>> resample_grid = np.array([1, 5, 9, 13, 14, 17, 21, 22, 23])
    >>> fluxc_resample = SplineInterpolatedResampler()
    >>> output_spectrum1D = fluxc_resample(input_spectra, resample_grid) # doctest: +IGNORE_OUTPUT

    """
    def __init__(self, bin_edges='nan_fill'):
        super().__init__(bin_edges)

    def __call__(self, orig_spectrum, fin_lamb):
        """
        Return the resulting `~specutils.Spectrum1D` of the resampling.
        """
        return self.resample1d(orig_spectrum, fin_lamb)

    def _interpolation(self, orig_dispersion, flux, fin_lamb):
        """
        Use specified interpolation to calculated resampled
        flux.

        Parameters
        ----------
        orig_dispersion : ndarray
            The original dispersion array.
        flux: ndarray
            The flux array from the input Spectrum1D
        fin_lamb : ndarray
            The desired dispersion array.

        Returns
        -------
        resample_flux : ndarray
            The resampled flux array generated from the interpolation.

        """
        cubic_spline = CubicSpline(orig_dispersion, flux, extrapolate=False)
        return cubic_spline(fin_lamb)

    def resample1d(self, orig_spectrum, fin_lamb):
        """
        Call interpolation, repackage new spectra


        Parameters
        ----------
        orig_spectrum : `~specutils.Spectrum1D`
            The original 1D spectrum.
        fin_lamb : ndarray
            The desired dispersion array.

        Returns
        -------
        resample_spectrum : `~specutils.Spectrum1D`
            An output spectrum containing the resampled `~specutils.Spectrum1D`
        """
        out_flux = self._interpolation(orig_spectrum.spectral_axis, orig_spectrum.flux,
                                  fin_lamb)

        # todo: for now, use the units from the pre-resampled
        # spectra, although if a unit is defined for fin_lamb and it doesn't
        # match the input spectrum it won't work right, will have to think
        # more about how to handle that... could convert before and after
        # calculation, which is probably easiest. Matrix math algorithm is
        # geometry based, so won't work to just let quantity math handle it.
        # todo: handle uncertainties for interpolated cases.
        resampled_spectrum = Spectrum1D(flux=out_flux * orig_spectrum.flux.unit,
                                        spectral_axis=np.array(fin_lamb) * orig_spectrum.spectral_axis_unit)

        return resampled_spectrum
