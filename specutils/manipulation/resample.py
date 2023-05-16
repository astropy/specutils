from abc import ABC, abstractmethod

import numpy as np
from astropy.nddata import VarianceUncertainty, InverseVariance
from astropy.units import Quantity
from scipy.interpolate import CubicSpline

from ..spectra import Spectrum1D, SpectralAxis

__all__ = ['ResamplerBase', 'FluxConservingResampler',
           'LinearInterpolatedResampler', 'SplineInterpolatedResampler']


class ResamplerBase(ABC):
    """
    Base class for resample classes.  The algorithms and needs for difference
    resamples will vary quite a bit, so this class is relatively sparse.

    Parameters
    ----------
    extrapolation_treatment : str
        What to do when resampling off the edge of the spectrum.  Can be
        ``'nan_fill'`` to have points beyond the edges by set to NaN, or
        ``'zero_fill'`` to be set to zero.
    """
    def __init__(self, extrapolation_treatment='nan_fill'):
        if extrapolation_treatment not in ('nan_fill', 'zero_fill'):
            raise ValueError('invalid extrapolation_treatment value: ' + str(extrapolation_treatment))
        self.extrapolation_treatment = extrapolation_treatment

    def __call__(self, orig_spectrum, fin_spec_axis):
        """
        Return the resulting `~specutils.Spectrum1D` of the resampling.
        """
        return self.resample1d(orig_spectrum, fin_spec_axis)

    @abstractmethod
    def resample1d(self, orig_spectrum, fin_spec_axis):
        """
        Workhorse method that will return the resampled Spectrum1D
        object.
        """
        return NotImplemented


class FluxConservingResampler(ResamplerBase):

    """
    This resampling algorithm conserves overall integrated flux (as opposed to
    flux density).
    Algorithm based on the equations documented in the following paper:
    https://ui.adsabs.harvard.edu/abs/2017arXiv170505165C/abstract

    Parameters
    ----------
    extrapolation_treatment : str
        What to do when resampling off the edge of the spectrum.  Can be
        ``'nan_fill'`` to have points beyond the edges by set to NaN, or
        ``'zero_fill'`` to be set to zero.

    Examples
    --------

    To resample an input spectrum to a user specified spectral grid using
    a flux conserving algorithm:

    >>> import numpy as np
    >>> import astropy.units as u
    >>> from specutils import Spectrum1D
    >>> from specutils.manipulation import FluxConservingResampler
    >>> input_spectra = Spectrum1D(
    ...     flux=np.array([1, 3, 7, 6, 20]) * u.mJy,
    ...     spectral_axis=np.array([2, 4, 12, 16, 20]) * u.nm)
    >>> resample_grid = [1, 5, 9, 13, 14, 17, 21, 22, 23]  *u.nm
    >>> fluxc_resample = FluxConservingResampler()
    >>> fluxc_resample(input_spectra, resample_grid)  # doctest: +FLOAT_CMP
    <Spectrum1D(flux=<Quantity [  nan,  3.  ,  6.  ,  7.  ,  6.25, 10.  , 20.  ,   nan,   nan] mJy>, spectral_axis=<SpectralAxis [ 1.,  5.,  9., 13., 14., 17., 21., 22., 23.] nm>)>
    """

    @staticmethod
    def _check_if_overlap(a1, a2, b1, b2):
        # check if 2 bins bounded by ``bin1_bounds``=(a1, a2) and
        # ``bin2_bounds``=(b1, b2) have any intersection
        if (b1 >= a2) or (b2 <= a1):
            return False
        return True

    @staticmethod
    def _get_overlap_area(a1, a2, b1, b2):
        # get intersection area of 2 bins bounded by
        # ``bin1_bounds``=(a1, b1) and ``bin2_bounds``=(a2, b2)
        # assuming its already known they overlap
        if b1 >= a1:
            if b2 > a2:
                return a2 - b1
            return b2 - b1
        if b2 > a2:
            return min(b2, a2) - max(b1, a1)
        return b2 - a1

    def _fluxc_resample(self, input_bin_centers, output_bin_centers,
                        input_bin_fluxes, errs):
        """
        Resample ``input_bin_fluxes`` and (optionally) ``errs`` from
        ``input_bin_centers`` to ``output_bin_centers``.

        Parameters
        ----------
        input_bin_centers : `~specutils.SpectralAxis`
            `~specutils.SpectralAxis` object, with input bin centers.
        output_bin_centers : `~specutils.SpectralAxis`
            `~specutils.SpectralAxis` object, with input bin centers.
        input_bin_fluxes : Quantity
            Quantity array of flux values.
        errs : `~astropy.nddata.Variance` object, or None
            Variance array of errors corresponding to input bin fluxes. If None,
            error resampling is not performed.

        Returns
        -------
       (output_fluxes, output_errs)
            A tuple containing plain numpy arrays of the resampled fluxes and
            errors (if available).
        """

        fill_val = np.nan  # bin_edges=nan_fill case
        if self.extrapolation_treatment == 'zero_fill':
            fill_val = 0

        # get bin edges from centers
        input_bin_edges = input_bin_centers.bin_edges.value
        output_bin_edges = output_bin_centers.bin_edges.value

        # convert to regular arrays
        input_bin_centers = input_bin_centers.value
        output_bin_centers = output_bin_centers.value
        input_bin_fluxes = input_bin_fluxes.value

        # create array of output fluxes and errors to be returned
        output_fluxes = np.zeros(shape=len(output_bin_centers))
        output_errs = None
        if errs is not None:
            output_errs = np.zeros(shape=len(output_bin_centers))

        # running sums for ouput flux, err calculations
        num = 0
        num_err = 0
        denom = 0

        # index of min. overlapping bin. - ignore bins before this in each iter
        exclude_before_last_iter = 0

        # first, figure out what output bins cover wavelengths outside the span of
        # input bins. these bins should have fluxes set to nan (or whatever the
        # fill val is.) and can be skipped
        min_idx = 0
        max_idx = None

        low_out_of_range = np.where(output_bin_edges < input_bin_edges[0])[0]
        if len(low_out_of_range) > 0:  # if any bins below wavelength range
            min_idx = low_out_of_range[-1] + 1
            output_fluxes[:min_idx] = fill_val
            if errs is not None:
                output_errs[:min_idx] = fill_val

        high_out_of_range = np.where(output_bin_edges > input_bin_edges[-1])[0]
        if len(high_out_of_range) > 0:
            max_idx = high_out_of_range[0] - 1
            output_fluxes[max_idx:] = fill_val
            if errs is not None:
                output_errs[max_idx:] = fill_val

        # iterate over each output bin in wavelength range of input bins
        for i, output_bin in enumerate(output_bin_centers[min_idx:max_idx]):

            i = i + min_idx  # index in orig, unclipped array
            bin_start, bin_stop = output_bin_edges[i], output_bin_edges[i+1]

            # now iterate over each input bin to see which ones overlap
            # but do this in such a way that we don't check bins unnecessarily.
            exclude_before = exclude_before_last_iter
            for j, input_bin in enumerate(input_bin_centers[exclude_before_last_iter:]):

                j = j + exclude_before  # index in orig, unclipped

                input_bin_start, input_bin_stop = input_bin_edges[j], input_bin_edges[j+1]
                input_bin_width = input_bin_stop - input_bin_start

                if input_bin_start >= bin_stop:  # stop when bins go out of range
                    break

                # if the start or the stop of the bin is in range, it overlaps
                overlaps = self._check_if_overlap(*(bin_start, bin_stop),
                                                  *(input_bin_start, input_bin_stop))

                if overlaps is True:

                    # if this input bin overlaps and terminates within this output
                    # bin, then it can be skipped for the next output bin (i)
                    if input_bin_stop <= bin_stop:
                        exclude_before_last_iter = j

                    # what fraction of input bin is contained in output bin?
                    overlap_width = self._get_overlap_area(*(bin_start, bin_stop),
                                                           *(input_bin_start,
                                                             input_bin_stop))
                    overlap_fraction = (overlap_width / input_bin_width)

                    # numerator contribution from this i, j
                    num += (input_bin_fluxes[j] * overlap_fraction * input_bin_width)
                    if errs is not None:
                        num_err += (overlap_fraction*overlap_fraction) *\
                                   (input_bin_width*input_bin_width) *\
                                   (errs[j]*errs[j])

                    # denom. contribution from this i, j
                    denom += input_bin_width * overlap_fraction

                else:  # bin doesn't overlap
                    continue

            output_fluxes[i] = num / denom

            if errs is not None:
                output_errs[i] = num_err / (denom * denom)
                num_err = 0  # reset running sum of P_ij^2*w_i^2*err_i^2

            num = 0  # reset running sum of P_ij*w_i*f_i
            denom = 0  # reset running sum of P_ij*w_i

        if errs is not None:
            output_errs = InverseVariance(np.reciprocal(np.sqrt(output_errs)))

        return (output_fluxes, output_errs)

    def resample1d(self, orig_spectrum, fin_spec_axis):
        """
        Create a re-sampling matrix to be used in re-sampling spectra in a way
        that conserves flux. If an uncertainty is present in the input spectra
        it will be propagated through to the final resampled output spectra
        as an InverseVariance uncertainty.

        Parameters
        ----------
        orig_spectrum : `~specutils.Spectrum1D`
            The original 1D spectrum.
        fin_spec_axis :  Quantity
            The desired spectral axis array.

        Returns
        -------
        resampled_spectrum : `~specutils.Spectrum1D`
            An output spectrum containing the resampled `~specutils.Spectrum1D`
        """

        if isinstance(fin_spec_axis, Quantity):
            if orig_spectrum.spectral_axis.unit != fin_spec_axis.unit:
                raise ValueError("Original spectrum spectral axis grid and new"
                                 "spectral axis grid must have the same units.")

        if not isinstance(fin_spec_axis, SpectralAxis):
            fin_spec_axis = SpectralAxis(fin_spec_axis)

        # Get provided uncertainty into variance
        if orig_spectrum.uncertainty is not None:
            pixel_uncer = orig_spectrum.uncertainty.represent_as(VarianceUncertainty).array
        else:
            pixel_uncer = None

        # convert unit
        orig_axis_in_fin = orig_spectrum.spectral_axis.to(fin_spec_axis.unit)

        # handle 2d flux inputs
        if orig_spectrum.flux.ndim == 2:
            # make output matrix
            output_fluxes = np.zeros(shape=(orig_spectrum.flux.shape[0],
                                            fin_spec_axis.shape[0]))
            output_errs = np.zeros(shape=(orig_spectrum.flux.shape[0],
                                          fin_spec_axis.shape[0]))
            for r, row in enumerate(orig_spectrum.flux):
                new_f, new_e = self._fluxc_resample(input_bin_centers=orig_axis_in_fin,
                                                    output_bin_centers=fin_spec_axis,
                                                    input_bin_fluxes=row,
                                                    errs=pixel_uncer[r])
                output_fluxes[r] = new_f
                output_errs[r] = new_e.array

            new_errs = InverseVariance(output_errs)

        else:
            # calculate new fluxes and errors
            output_fluxes, new_errs = self._fluxc_resample(input_bin_centers=orig_axis_in_fin,
                                                           output_bin_centers=fin_spec_axis,
                                                           input_bin_fluxes=orig_spectrum.flux,
                                                           errs=pixel_uncer)

        output_fluxes = output_fluxes << orig_spectrum.flux.unit
        fin_spec_axis = np.array(fin_spec_axis) << orig_spectrum.spectral_axis.unit

        resampled_spectrum = Spectrum1D(flux=output_fluxes,
                                        spectral_axis=fin_spec_axis,
                                        uncertainty=new_errs)

        return resampled_spectrum


class LinearInterpolatedResampler(ResamplerBase):
    """
    Resample a spectrum onto a new ``spectral_axis`` using linear interpolation.

    Parameters
    ----------
    extrapolation_treatment : str
        What to do when resampling off the edge of the spectrum.  Can be
        ``'nan_fill'`` to have points beyond the edges by set to NaN, or
        ``'zero_fill'`` to be set to zero.

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
    >>> resample_grid = [1, 5, 9, 13, 14, 17, 21, 22, 23] * u.nm
    >>> fluxc_resample = LinearInterpolatedResampler()
    >>> fluxc_resample(input_spectra, resample_grid)  # doctest: +FLOAT_CMP
    <Spectrum1D(flux=<Quantity [ nan, 3.5 , 5.5 , 6.75, 6.5 , 9.5 ,  nan,  nan,  nan] mJy>, spectral_axis=<SpectralAxis [ 1.,  5.,  9., 13., 14., 17., 21., 22., 23.] nm>)>

    """
    def __init__(self, extrapolation_treatment='nan_fill'):
        super().__init__(extrapolation_treatment)

    def resample1d(self, orig_spectrum, fin_spec_axis):
        """
        Call interpolation, repackage new spectra


        Parameters
        ----------
        orig_spectrum : `~specutils.Spectrum1D`
            The original 1D spectrum.
        fin_spec_axis : ndarray
            The desired spectral axis array.

        Returns
        -------
        resample_spectrum : `~specutils.Spectrum1D`
            An output spectrum containing the resampled `~specutils.Spectrum1D`
        """

        fill_val = np.nan  # bin_edges=nan_fill case
        if self.extrapolation_treatment == 'zero_fill':
            fill_val = 0

        orig_axis_in_fin = orig_spectrum.spectral_axis.to(fin_spec_axis.unit)

        out_flux_arr = np.interp(fin_spec_axis.value, orig_axis_in_fin.value,
                                 orig_spectrum.flux.value, left=fill_val, right=fill_val)
        out_flux = Quantity(out_flux_arr, unit=orig_spectrum.flux.unit)

        new_unc = None
        if orig_spectrum.uncertainty is not None:
            out_unc_arr = np.interp(fin_spec_axis.value, orig_axis_in_fin.value,
                                    orig_spectrum.uncertainty.array,
                                    left=fill_val, right=fill_val)
            new_unc = orig_spectrum.uncertainty.__class__(array=out_unc_arr,
                                                          unit=orig_spectrum.unit)

        return Spectrum1D(spectral_axis=fin_spec_axis,
                          flux=out_flux,
                          uncertainty=new_unc)


class SplineInterpolatedResampler(ResamplerBase):
    """
    This resample algorithim uses a cubic spline interpolator. Any uncertainty
    is also interpolated using an identical spline.


    Parameters
    ----------
    extrapolation_treatment : str
        What to do when resampling off the edge of the spectrum.  Can be
        ``'nan_fill'`` to have points beyond the edges by set to NaN, or
        ``'zero_fill'`` to be set to zero.

    Examples
    --------

    To resample an input spectrum to a user specified spectral axis grid using
    a cubic spline interpolator:

    >>> import numpy as np
    >>> import astropy.units as u
    >>> from specutils import Spectrum1D
    >>> from specutils.manipulation import SplineInterpolatedResampler
    >>> input_spectra = Spectrum1D(
    ...     flux=np.array([1, 3, 7, 6, 20]) * u.mJy,
    ...     spectral_axis=np.array([2, 4, 12, 16, 20]) * u.nm)
    >>> resample_grid = [1, 5, 9, 13, 14, 17, 21, 22, 23] * u.nm
    >>> fluxc_resample = SplineInterpolatedResampler()
    >>> fluxc_resample(input_spectra, resample_grid)  # doctest: +FLOAT_CMP
    <Spectrum1D(flux=<Quantity [       nan, 3.98808594, 6.94042969, 6.45869141, 5.89921875,
               7.29736328,        nan,        nan,        nan] mJy>, spectral_axis=<SpectralAxis [ 1.,  5.,  9., 13., 14., 17., 21., 22., 23.] nm>)>

    """
    def __init__(self, bin_edges='nan_fill'):
        super().__init__(bin_edges)

    def resample1d(self, orig_spectrum, fin_spec_axis):
        """
        Call interpolation, repackage new spectra


        Parameters
        ----------
        orig_spectrum : `~specutils.Spectrum1D`
            The original 1D spectrum.
        fin_spec_axis : Quantity
            The desired spectral axis array.

        Returns
        -------
        resample_spectrum : `~specutils.Spectrum1D`
            An output spectrum containing the resampled `~specutils.Spectrum1D`
        """
        orig_axis_in_new = orig_spectrum.spectral_axis.to(fin_spec_axis.unit)
        flux_spline = CubicSpline(orig_axis_in_new.value, orig_spectrum.flux.value,
                                  extrapolate=self.extrapolation_treatment != 'nan_fill')
        out_flux_val = flux_spline(fin_spec_axis.value)

        new_unc = None
        if orig_spectrum.uncertainty is not None:
            unc_spline = CubicSpline(orig_axis_in_new.value, orig_spectrum.uncertainty.array,
                                     extrapolate=self.extrapolation_treatment != 'nan_fill')
            out_unc_val = unc_spline(fin_spec_axis.value)
            new_unc = orig_spectrum.uncertainty.__class__(array=out_unc_val,
                                                          unit=orig_spectrum.unit)

        if self.extrapolation_treatment == 'zero_fill':
            origedges = orig_spectrum.spectral_axis.bin_edges
            off_edges = (fin_spec_axis < origedges[0]) | (origedges[-1] < fin_spec_axis)
            out_flux_val[off_edges] = 0
            if new_unc is not None:
                new_unc.array[off_edges] = 0

        return Spectrum1D(spectral_axis=fin_spec_axis,
                          flux=out_flux_val*orig_spectrum.flux.unit,
                          uncertainty=new_unc)
