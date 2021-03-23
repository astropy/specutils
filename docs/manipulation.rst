====================
Manipulating Spectra
====================

While there are myriad ways you might want to alter a spectrum,
:ref:`specutils <specutils>` provides some specific functionality that
is commonly used in astronomy.  These tools are detailed here, but it
is important to bear in mind that this is *not* intended to be
exhaustive - the point of :ref:`specutils <specutils>` is to provide a
framework you can use to do your data analysis.  Hence the
functionality described here is best thought of as pieces you might
string together with your own functionality to build a tailor-made
spectral analysis environment.

In general, however, :ref:`specutils <specutils>` is designed around
the idea that spectral manipulations generally yield *new* spectrum
objects, rather than in-place operations.  This is not a true
restriction, but is a guideline that is recommended primarily to keep
you from accidentally modifying a spectrum you didn't mean to change.

Smoothing
---------

Specutils provides smoothing for spectra in two forms: 1) convolution based
using smoothing `astropy.convolution` and 2) median filtering
using the :func:`scipy.signal.medfilt`.  Each of these act on the flux
of the :class:`~specutils.Spectrum1D` object.

.. note:: Specutils smoothing kernel widths and standard deviations are
             in units of pixels and not ``Quantity``.

Convolution Based Smoothing
^^^^^^^^^^^^^^^^^^^^^^^^^^^

While any kernel supported by `astropy.convolution` will work (using the
`~specutils.manipulation.convolution_smooth` function), several
commonly-used kernels have convenience functions wrapping them to simplify
the smoothing process into a simple one-line operation.  Currently
implemented are: :func:`~specutils.manipulation.box_smooth`
(:class:`~astropy.convolution.Box1DKernel`),
:func:`~specutils.manipulation.gaussian_smooth`
(:class:`~astropy.convolution.Gaussian1DKernel`), and
:func:`~specutils.manipulation.trapezoid_smooth`
(:class:`~astropy.convolution.Trapezoid1DKernel`).


.. code-block:: python

    >>> from specutils import Spectrum1D
    >>> import astropy.units as u
    >>> import numpy as np
    >>> from specutils.manipulation import (box_smooth, gaussian_smooth, trapezoid_smooth)

    >>> spec1 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm, flux=np.random.sample(49)*u.Jy)
    >>> spec1_bsmooth = box_smooth(spec1, width=3)
    >>> spec1_gsmooth = gaussian_smooth(spec1, stddev=3)
    >>> spec1_tsmooth = trapezoid_smooth(spec1, width=3)
    >>> gaussian_smooth(spec1, stddev=3) #doctest:+SKIP
    Spectrum1D([0.22830748, 0.2783204 , 0.32007408, 0.35270403, 0.37899655,
                0.40347983, 0.42974259, 0.45873436, 0.48875214, 0.51675647,
                0.54007149, 0.55764758, 0.57052796, 0.58157173, 0.59448669,
                0.61237409, 0.63635755, 0.66494062, 0.69436655, 0.7199299 ,
                0.73754271, 0.74463192, 0.74067744, 0.72689092, 0.70569365,
                0.6800534 , 0.65262146, 0.62504013, 0.59778884, 0.57072578,
                0.54416776, 0.51984003, 0.50066938, 0.48944714, 0.48702192,
                0.49126444, 0.49789092, 0.50276877, 0.50438924, 0.50458914,
                0.50684731, 0.51321106, 0.52197328, 0.52782086, 0.52392599,
                0.50453064, 0.46677128, 0.41125485, 0.34213489])


Each of the specific smoothing methods create the appropriate `astropy.convolution.convolve`
kernel and then call a helper function :func:`~specutils.manipulation.convolution_smooth`
that takes the spectrum and an astropy 1D kernel.  So, one could also do:

.. code-block:: python

    >>> from astropy.convolution import Box1DKernel
    >>> from specutils.manipulation import convolution_smooth

    >>> box1d_kernel = Box1DKernel(width=3)

    >>> spec1 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm, flux=np.random.sample(49) * u.Jy)
    >>> spec1_bsmooth2 = convolution_smooth(spec1, box1d_kernel)

In this case, the ``spec1_bsmooth2`` result should be equivalent to the ``spec1_bsmooth`` in
the section above (assuming the flux data of the input ``spec`` is the same).

The uncertainties are propagated using a standard "propagation of errors" method, if the uncertainty
is defined for the spectrum *and* it is one of StdDevUncertainty, VarianceUncertainty or InverseVariance.
But note that this does *not* consider covariance between points.

Median Smoothing
^^^^^^^^^^^^^^^^

The median based smoothing  is implemented using `scipy.signal.medfilt` and
has a similar call structure to the convolution-based smoothing methods. This
method applys the median filter across the flux.


.. note::
    This method is not flux conserving and errors are not propagated.


.. code-block:: python

    >>> from specutils.manipulation import median_smooth

    >>> spec1 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm, flux=np.random.sample(49) * u.Jy)
    >>> spec1_msmooth = median_smooth(spec1, width=3)

Resampling
----------
:ref:`specutils <specutils>` contains several classes for resampling the flux
in a :class:`~specutils.Spectrum1D` object.  Currently supported methods of
resampling are integrated flux conserving with :class:`~specutils.manipulation.FluxConservingResampler`,
linear interpolation with :class:`~specutils.manipulation.LinearInterpolatedResampler`,
and cubic spline with :class:`~specutils.manipulation.SplineInterpolatedResampler`.
Each of these classes takes in a :class:`~specutils.Spectrum1D` and a user
defined output dispersion grid, and returns a new :class:`~specutils.Spectrum1D`
with the resampled flux. Currently the resampling classes expect the new
dispersion grid unit to be the same as the input spectrum's dispersion grid unit.

If the input :class:`~specutils.Spectrum1D` contains an uncertainty,
:class:`~specutils.manipulation.FluxConservingResampler` will propogate the
uncertainty to the final output :class:`~specutils.Spectrum1D`. However, the
other two implemented resampling classes (:class:`~specutils.manipulation.LinearInterpolatedResampler`
and :class:`~specutils.manipulation.SplineInterpolatedResampler`) will ignore
any input uncertainty.

Here's a set of simple examples showing each of the three types of resampling:

.. plot::
    :include-source:
    :align: center
    :context: close-figs

    First are the imports we will need as well as loading in the example data:

    >>> from astropy.io import fits
    >>> from astropy import units as u
    >>> import numpy as np
    >>> from matplotlib import pyplot as plt
    >>> from astropy.visualization import quantity_support
    >>> quantity_support()  # for getting units on the axes below  # doctest: +IGNORE_OUTPUT

    >>> f = fits.open('https://data.sdss.org/sas/dr16/sdss/spectro/redux/26/spectra/1323/spec-1323-52797-0012.fits')  # doctest: +IGNORE_OUTPUT +REMOTE_DATA
    >>> # The spectrum is in the second HDU of this file.
    >>> specdata = f[1].data[1020:1250] # doctest: +REMOTE_DATA
    >>> f.close() # doctest: +REMOTE_DATA

    Then we re-format this dataset into astropy quantities, and create a
    `~specutils.Spectrum1D` object:

    >>> from specutils import Spectrum1D
    >>> lamb = 10**specdata['loglam'] * u.AA # doctest: +REMOTE_DATA
    >>> flux = specdata['flux'] * 10**-17 * u.Unit('erg cm-2 s-1 AA-1') # doctest: +REMOTE_DATA
    >>> input_spec = Spectrum1D(spectral_axis=lamb, flux=flux) # doctest: +REMOTE_DATA

    >>> f, ax = plt.subplots()  # doctest: +IGNORE_OUTPUT
    >>> ax.step(input_spec.spectral_axis, input_spec.flux) # doctest: +IGNORE_OUTPUT +REMOTE_DATA

.. plot::
    :include-source:
    :align: center
    :context: close-figs

    Now we show examples and plots of the different resampling currently
    available.

    >>> from specutils.manipulation import FluxConservingResampler, LinearInterpolatedResampler, SplineInterpolatedResampler
    >>> new_disp_grid = np.arange(4800, 5200, 3) * u.AA

    Flux Conserving Resampler:

    >>> fluxcon = FluxConservingResampler()
    >>> new_spec_fluxcon = fluxcon(input_spec, new_disp_grid) # doctest: +IGNORE_OUTPUT +REMOTE_DATA
    >>> f, ax = plt.subplots()  # doctest: +IGNORE_OUTPUT
    >>> ax.step(new_spec_fluxcon.spectral_axis, new_spec_fluxcon.flux) # doctest: +IGNORE_OUTPUT +REMOTE_DATA

.. plot::
    :include-source:
    :align: center
    :context: close-figs

    Linear Interpolation Resampler:

    >>> linear = LinearInterpolatedResampler()
    >>> new_spec_lin = linear(input_spec, new_disp_grid)  # doctest: +REMOTE_DATA
    >>> f, ax = plt.subplots()  # doctest: +IGNORE_OUTPUT
    >>> ax.step(new_spec_lin.spectral_axis, new_spec_lin.flux) # doctest: +IGNORE_OUTPUT +REMOTE_DATA

.. plot::
    :include-source:
    :align: center
    :context: close-figs

    Spline Resampler:

    >>> spline = SplineInterpolatedResampler()
    >>> new_spec_sp = spline(input_spec, new_disp_grid)  # doctest: +REMOTE_DATA
    >>> f, ax = plt.subplots()  # doctest: +IGNORE_OUTPUT
    >>> ax.step(new_spec_sp.spectral_axis, new_spec_sp.flux) # doctest: +IGNORE_OUTPUT +REMOTE_DATA

Splicing/Combining Multiple Spectra
-----------------------------------
The resampling functionality detailed above is also the default way
:ref:`specutils <specutils>` supports splicing multiple spectra together into a
single spectrum. This can be achieved as follows:

.. plot::
    :include-source:
    :align: center
    :context: close-figs

    >>> spec1 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.micron, flux=np.random.randn(49)*u.Jy)
    >>> spec2 = Spectrum1D(spectral_axis=np.arange(51, 100) * u.micron, flux=(np.random.randn(49)+1)*u.Jy)

    >>> new_spectral_axis = np.concatenate([spec1.spectral_axis.value, spec2.spectral_axis.to_value(spec1.spectral_axis.unit)]) * spec1.spectral_axis.unit

    >>> resampler = LinearInterpolatedResampler(extrapolation_treatment='zero_fill')
    >>> new_spec1 = resampler(spec1, new_spectral_axis)
    >>> new_spec2 = resampler(spec2, new_spectral_axis)

    >>> final_spec = new_spec1 + new_spec2

    Yielding a spliced spectrum (the solid line below) composed of the splice of
    two other spectra (dashed lines)::

    >>> f, ax = plt.subplots()  # doctest: +IGNORE_OUTPUT
    >>> ax.step(final_spec.spectral_axis, final_spec.flux, where='mid', c='k', lw=2) # doctest: +IGNORE_OUTPUT
    >>> ax.step(spec1.spectral_axis, spec1.flux, ls='--', where='mid', lw=1) # doctest: +IGNORE_OUTPUT
    >>> ax.step(spec2.spectral_axis, spec2.flux, ls='--', where='mid', lw=1) # doctest: +IGNORE_OUTPUT



Uncertainty Estimation
----------------------

Some of the machinery in :ref:`specutils <specutils>` (e.g.
`~specutils.analysis.snr`) requires an uncertainty to be present.
While some data reduction pipelines generate this as part of the
reduction process, sometimes it's necessary to estimate the
uncertainty in a spectrum using the spectral data itself. Currently
:ref:`specutils <specutils>` provides the straightforward
`~specutils.manipulation.noise_region_uncertainty` function.

First we build a spectrum like that used in :doc:`analysis`, but without a
known uncertainty:

.. code-block:: python

    >>> from astropy.modeling import models
    >>> np.random.seed(42)
    >>> spectral_axis = np.linspace(10, 1, 200) * u.GHz
    >>> spectral_model = models.Gaussian1D(amplitude=3*u.Jy, mean=5*u.GHz, stddev=0.8*u.GHz)
    >>> flux = spectral_model(spectral_axis)
    >>> flux += np.random.normal(0., 0.2, spectral_axis.shape) * u.Jy
    >>> noisy_gaussian = Spectrum1D(spectral_axis=spectral_axis, flux=flux)

Now we estimate the uncertainty from the region that does *not* contain
the line:

.. code-block:: python

    >>> from specutils import SpectralRegion
    >>> from specutils.manipulation import noise_region_uncertainty
    >>> noise_region = SpectralRegion([(10, 7), (3, 0)] * u.GHz)
    >>> spec_w_unc = noise_region_uncertainty(noisy_gaussian, noise_region)
    >>> spec_w_unc.uncertainty # doctest: +ELLIPSIS
    StdDevUncertainty([0.1886835, ..., 0.1886835])

Or similarly, expressed in pixels:

.. code-block:: python

    >>> noise_region = SpectralRegion([(0, 25), (175, 200)]*u.pix)
    >>> spec_w_unc = noise_region_uncertainty(noisy_gaussian, noise_region)
    >>> spec_w_unc.uncertainty # doctest: +ELLIPSIS
    StdDevUncertainty([0.18739524, ..., 0.18739524])

S/N Threshold Mask
------------------

It is useful to be able to find all the spaxels in an ND spectrum
in which the signal to noise ratio is greater than some threshold.
This method implements this functionality so that a `~specutils.Spectrum1D`
object, `~specutils.SpectrumCollection` or an :class:`~astropy.nddata.NDData` derived
object may be passed in as the first parameter. The second parameter
is a floating point threshold.

For example, first a spectrum with flux and uncertainty is created, and
then call the ``snr_threshold`` method:

.. code-block:: python

    >>> import numpy as np
    >>> from astropy.nddata import StdDevUncertainty
    >>> import astropy.units as u
    >>> from specutils import Spectrum1D
    >>> from specutils.manipulation import snr_threshold
    >>> np.random.seed(42)
    >>> wavelengths = np.arange(0, 10)*u.um
    >>> flux = 100*np.abs(np.random.randn(10))*u.Jy
    >>> uncertainty = StdDevUncertainty(np.abs(np.random.randn(10))*u.Jy)
    >>> spectrum = Spectrum1D(spectral_axis=wavelengths, flux=flux, uncertainty=uncertainty)
    >>> spectrum_masked = snr_threshold(spectrum, 50) #doctest:+SKIP
    >>> # To create a masked flux array
    >>> flux_masked = spectrum_masked.flux #doctest:+SKIP
    >>> flux_masked[spectrum_masked.mask] = np.nan #doctest:+SKIP

The output ``spectrum_masked`` is a shallow copy of the input ``spectrum``
with the ``mask`` attribute set to False where the S/N is greater than 50
and True elsewhere. It is this way to be consistent with ``astropy.nddata``.

.. note:: The mask attribute is the only attribute modified by ``snr_threshold()``. To
             retrieve the masked flux data use ``spectrum.masked.flux_masked``.

Shifting
--------

In addition to resampling, you may sometimes wish to simply shift the
``spectral_axis`` of a spectrum (a la the ``specshift`` iraf task).
There is no explicit function for this because it is a basic transform of
the ``spectral_axis``. Therefore one can use a construct like this:

.. code-block:: python

    >>> from specutils import Spectrum1D
    >>> np.random.seed(42)
    >>> wavelengths = np.arange(0, 10) * u.um
    >>> flux = 100 * np.abs(np.random.randn(10)) * u.Jy
    >>> spectrum = Spectrum1D(spectral_axis=wavelengths, flux=flux)
    >>> spectrum  #doctest:+ELLIPSIS
    <Spectrum1D(flux=<Quantity [ 49.6714153 ,  13.82643012,  64.76885381, 152.30298564,
                23.41533747,  23.41369569, 157.92128155,  76.74347292,
                46.94743859,  54.25600436] Jy>, spectral_axis=<SpectralAxis [0., 1., 2., 3., 4., 5., 6., 7., 8., 9.] um>)>

    >>> shift = 12300 * u.AA
    >>> new_spec = Spectrum1D(spectral_axis=spectrum.spectral_axis + shift, flux=spectrum.flux)
    >>> new_spec #doctest:+ELLIPSIS
    <Spectrum1D(flux=<Quantity [ 49.6714153 ,  13.82643012,  64.76885381, 152.30298564,
                23.41533747,  23.41369569, 157.92128155,  76.74347292,
                46.94743859,  54.25600436] Jy>, spectral_axis=<SpectralAxis [ 1.23,  2.23,  3.23,  4.23,  5.23,  6.23,  7.23,  8.23,  9.23, 10.23] um>)>


Reference/API
-------------

.. automodapi:: specutils.manipulation
    :no-heading:
