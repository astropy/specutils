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
        ``'nan_fill'`` to have points beyond the edges by set to NaN,
        ``'zero_fill'`` to set thoe points to zero, or ``'truncate'`` to
        truncate any non-overlapping bins of the spectrum.
    """
    def __init__(self, extrapolation_treatment='nan_fill'):
        if extrapolation_treatment not in ('nan_fill', 'zero_fill', 'truncate'):
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
        ``'nan_fill'`` to have points beyond the edges by set to NaN,
        ``'zero_fill'`` to set those points to zero, or ``'truncate'`` to
        truncate any non-overlapping bins of the spectrum.

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
    <Spectrum1D(flux=<Quantity [ 1.  ,  3.  ,  6.  ,  7.  ,  6.25, 10.  , 20.  ,   nan,   nan] mJy> (shape=(9,), mean=7.60714 mJy); spectral_axis=<SpectralAxis [ 1.  5.  9. ... 21. 22. 23.] nm> (length=9))>

    """

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

        # create array of output fluxes and errors to be returned
        output_fluxes = np.zeros(shape=len(output_bin_centers)) * input_bin_fluxes.unit
        output_errs = None
        if errs is not None:
            output_errs = np.zeros(shape=len(output_bin_centers))

        # first, figure out what output bins cover wavelengths outside the span of
        # input bins. these bins should have fluxes set to nan (or whatever the
        # fill val is.) and can be skipped
        min_idx = 0
        max_idx = None

        low_out_of_range = np.where(output_bin_edges <= input_bin_edges[0])[0]
        if len(low_out_of_range) > 0:  # if any bins below wavelength range
            min_idx = low_out_of_range[-1]  # This doesn't need +1 because bin_edges has len+1 compared to output_fluxes
            output_fluxes[:min_idx] = fill_val
            if errs is not None:
                output_errs[:min_idx] = fill_val

        high_out_of_range = np.where(output_bin_edges > input_bin_edges[-1])[0]
        if len(high_out_of_range) > 0:
            max_idx = high_out_of_range[0] - 1
            output_fluxes[max_idx:] = fill_val
            if errs is not None:
                output_errs[max_idx:] = fill_val

        clipped_output_centers = output_bin_centers[min_idx:max_idx]

        # find the index of the first input bin that intersects the first
        # in-range output bin.
        first_output_edge = output_bin_edges[min_idx]

        idx_last_overlapping_bin = np.where(input_bin_edges[1:] > first_output_edge)[0][0]

        # iterate over each output bin in wavelength range of input bins
        for i, output_bin in enumerate(clipped_output_centers):

            i = i + min_idx  # index in orig, unclipped array

            bin_start, bin_stop = output_bin_edges[i], output_bin_edges[i+1]

            # the first at least partially overlapping bin was determined in the
            # last iteration (or by the initial clipping, if i=0)
            first_bin = idx_last_overlapping_bin

            # keep checking bins, starting at the one after we know overlaps first,
            # and stop when the back edge of an input bin overlaps the front
            # edge of this output bin.
            while input_bin_edges[idx_last_overlapping_bin + 1] < bin_stop:
                idx_last_overlapping_bin += 1

            # if the front edge of the last overlapping bin terminates in this
            # output bin, don't check it next time
            final_bin = idx_last_overlapping_bin

            if input_bin_edges[idx_last_overlapping_bin + 1] <= bin_stop:
                idx_last_overlapping_bin = idx_last_overlapping_bin + 1

            # now, calculate fluxes and errors

            # if only one input bin covers this output bin
            # flux_j=p_ij*w_i*f_i/p_ij*w_i = f_i - f_j=f_i, err_j=err_i
            if final_bin == first_bin:
                output_fluxes[i] = input_bin_fluxes[first_bin]
                if errs is not None:
                    output_errs[i] = errs[first_bin]

            # otherwise, figure out the contribution from each overlapping
            # input bin to calculate the final flux in the output bin.
            else:
                # the first edges of each overlapping input bin
                first_edges = input_bin_edges[first_bin:final_bin+1]
                # the final edges of each overlapping input bin
                second_edges = input_bin_edges[first_bin+1:final_bin+2]

                # to calculate the overlap area, of input on output, we
                # want to only deal with input bin's leading edges if they are
                # inside the output bin. otherwise, they are bounded by the
                # output bin edges. temporarily set the last edges to the output
                # bin bounds and then reset them at the end
                first_edges_orig_first = first_edges[0]
                first_edges[0] = bin_start
                second_edges_orig_last = second_edges[-1]
                second_edges[-1] = bin_stop

                p_ij = second_edges - first_edges

                # reset back
                first_edges[0] = first_edges_orig_first
                second_edges[-1] = second_edges_orig_last

                sum_pij = np.sum(p_ij)

                final_flux = (np.sum(input_bin_fluxes[first_bin:final_bin+1] * p_ij)) / sum_pij
                output_fluxes[i] = final_flux

                if errs is not None:
                    final_err = np.sum((errs[first_bin:final_bin+1] * p_ij) ** 2) / (sum_pij * sum_pij)
                    output_errs[i] = np.sqrt(final_err)

        if errs is not None:
            output_errs = InverseVariance(np.reciprocal(output_errs))

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

        # handle multi dimensional flux inputs
        if orig_spectrum.flux.ndim >= 2:

            # the output fluxes and errs should have the same shape as the input
            # except for the last axis, which should be the size of the new
            # spectral axis
            new_shape = tuple(list(orig_spectrum.shape[0:-1]) +
                              [len(fin_spec_axis)])

            # make output matricies
            output_fluxes = np.zeros(shape=new_shape)
            output_errs = np.zeros(shape=new_shape)

            for index, row in np.ndenumerate(orig_spectrum.flux[..., 0]):

                orig_fluxes = orig_spectrum.flux[index]
                orig_uncer = pixel_uncer[index]

                new_f, new_e = self._fluxc_resample(input_bin_centers=orig_axis_in_fin,
                                                    output_bin_centers=fin_spec_axis,
                                                    input_bin_fluxes=orig_fluxes,
                                                    errs=orig_uncer)
                output_fluxes[index] = new_f
                output_errs[index] = new_e.array

            new_errs = InverseVariance(output_errs)

        else:
            # calculate new fluxes and errors
            output_fluxes, new_errs = self._fluxc_resample(input_bin_centers=orig_axis_in_fin,
                                                           output_bin_centers=fin_spec_axis,
                                                           input_bin_fluxes=orig_spectrum.flux,
                                                           errs=pixel_uncer)

        output_fluxes = output_fluxes << orig_spectrum.flux.unit
        fin_spec_axis = np.array(fin_spec_axis) << orig_spectrum.spectral_axis.unit

        if self.extrapolation_treatment == 'truncate':
            fin_spec_axis = fin_spec_axis[np.where(~np.isnan(output_fluxes))]
            if new_errs is not None:
                new_errs = new_errs[np.where(~np.isnan(output_fluxes))]
            output_fluxes = output_fluxes[np.where(~np.isnan(output_fluxes))]

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
        ``'nan_fill'`` to have points beyond the edges by set to NaN,
        ``'zero_fill'`` to set those points to zero, or ``'truncate'`` to
        truncate any non-overlapping bins of the spectrum.

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
    <Spectrum1D(flux=<Quantity [ nan, 3.5 , 5.5 , 6.75, 6.5 , 9.5 ,  nan,  nan,  nan] mJy> (shape=(9,), mean=6.35000 mJy); spectral_axis=<SpectralAxis [ 1.  5.  9. ... 21. 22. 23.] nm> (length=9))>

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

        if self.extrapolation_treatment == 'truncate':
            fin_spec_axis = fin_spec_axis[np.where(~np.isnan(out_flux))]
            if new_unc is not None:
                new_unc = new_unc[np.where(~np.isnan(out_flux))]
            out_flux = out_flux[np.where(~np.isnan(out_flux))]

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
        ``'nan_fill'`` to have points beyond the edges by set to NaN,
        ``'zero_fill'`` to set those points to zero, or ``'truncate'`` to
        truncate any non-overlapping bins of the spectrum. Any other value will
        have the spline interpolate beyond the edges of the original data.

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
               7.29736328,        nan,        nan,        nan] mJy> (shape=(9,), mean=6.11676 mJy); spectral_axis=<SpectralAxis [ 1.  5.  9. ... 21. 22. 23.] nm> (length=9))>

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
        fin_spec_axis : Quantity
            The desired spectral axis array.

        Returns
        -------
        resample_spectrum : `~specutils.Spectrum1D`
            An output spectrum containing the resampled `~specutils.Spectrum1D`
        """
        orig_axis_in_new = orig_spectrum.spectral_axis.to(fin_spec_axis.unit)
        flux_spline = CubicSpline(orig_axis_in_new.value, orig_spectrum.flux.value,
                                  extrapolate=self.extrapolation_treatment not in ('nan_fill',
                                                                                   'zero_fill',
                                                                                   'truncate'))
        out_flux_val = flux_spline(fin_spec_axis.value)

        new_unc = None
        if orig_spectrum.uncertainty is not None:
            unc_spline = CubicSpline(orig_axis_in_new.value, orig_spectrum.uncertainty.array,
                                     extrapolate=self.extrapolation_treatment not in ('nan_fill',
                                                                                      'zero_fill',
                                                                                      'truncate'))
            out_unc_val = unc_spline(fin_spec_axis.value)
            new_unc = orig_spectrum.uncertainty.__class__(array=out_unc_val, unit=orig_spectrum.unit)

        fill_val = np.nan
        if self.extrapolation_treatment == 'zero_fill':
            fill_val = 0

        origedges = orig_spectrum.spectral_axis.bin_edges
        off_edges = (fin_spec_axis < origedges[0]) | (origedges[-1] < fin_spec_axis)
        out_flux_val[off_edges] = fill_val
        if new_unc is not None:
            new_unc.array[off_edges] = fill_val

        if self.extrapolation_treatment == 'truncate':
            fin_spec_axis = fin_spec_axis[np.where(~np.isnan(out_flux_val))]
            if new_unc is not None:
                new_unc = new_unc[np.where(~np.isnan(out_flux_val))]
            out_flux_val = out_flux_val[np.where(~np.isnan(out_flux_val))]

        return Spectrum1D(spectral_axis=fin_spec_axis,
                          flux=out_flux_val*orig_spectrum.flux.unit,
                          uncertainty=new_unc)
