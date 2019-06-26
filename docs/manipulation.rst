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


Median Smoothing
^^^^^^^^^^^^^^^^

The median based smoothing  is implemented using `scipy.signal.medfilt` and
has a similar call structure to the convolution-based smoothing methods. This
method applys the median filter across the flux.

Note: This method is not flux conserving.

.. code-block:: python

    >>> from specutils.manipulation import median_smooth

    >>> spec1 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm, flux=np.random.sample(49) * u.Jy)
    >>> spec1_msmooth = median_smooth(spec1, width=3)

Resampling
----------
:ref:`specutils <specutils>` contains several classes for resampling the flux
in a :class:`~specutils.Spectrum1D` object.  Currently supported methods of
resampling are integrated flux conserving with :class:`~specutils.manipulation.FluxConservingResample`,
linear interpolation with :class:`~specutils.manipulation.LinearInterpolatedResample`,
and cubic spline with :class:`~specutils.manipulation.SplineInterpolatedResample`.
Each of these classes takes in a :class:`~specutils.Spectrum1D` and a user
defined output dispersion grid, and returns a new :class:`~specutils.Spectrum1D`
with the resampled flux. Currently the resampling classes expect the new
dispersion grid unit to be the same as the input spectrum's dispersion grid unit.

If the input :class:`~specutils.Spectrum1D` contains an uncertainty,
:class:`~specutils.manipulation.FluxConservingResample` will propogate the
uncertainty to the final output :class:`~specutils.Spectrum1D`. However, the
other two implemented resampling classes (:class:`~specutils.manipulation.LinearInterpolatedResample`
and :class:`~specutils.manipulation.SplineInterpolatedResample`) will ignore
any input uncertainty.

Here's a set of simple examples showing each of the three types of resampling:

.. code-block:: python

    >>> from specutils import Spectrum1D
    >>> import astropy.units as u
    >>> import numpy as np
    >>> from specutils.manipulation import FluxConservingResample, LinearInterpolatedResample, SplineInterpolatedResample

    >>> input_spec = Spectrum1D(spectral_axis=np.arange(1, 20, 2) * u.nm, flux=np.random.sample(10)*u.Jy)
    >>> new_disp_grid = np.arange(1.5,20.5)

    >>> fluxcon = FluxConservingResample()
    >>> new_spec = fluxcon(input_spec, new_disp_grid)
    >>> new_spec.spectral_axis #doctest:+SKIP
    <Quantity [ 1.5,  2.5,  3.5,  4.5,  5.5,  6.5,  7.5,  8.5,  9.5, 10.5,
               11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5] nm>

    >>> new_spec.flux #doctest:+SKIP
    <Quantity [0.51712972, 0.10433299, 0.10433299, 0.8832977 , 0.8832977 ,
               0.9894764 , 0.9894764 , 0.25040808, 0.25040808, 0.10664375,
               0.10664375, 0.58550929, 0.58550929, 0.50498475, 0.50498475,
               0.42223018, 0.42223018, 0.75861711, 0.75861711] Jy>

    >>> linear = LinearInterpolatedResample()
    >>> new_spec2 = linear(input_spec, new_disp_grid)
    >>> new_spec2.flux #doctest:+SKIP
    <Quantity [0.41393054, 0.20753217, 0.29907417, 0.68855652, 0.90984238,
               0.96293173, 0.80470932, 0.43517516, 0.214467  , 0.14258483,
               0.22636014, 0.4657929 , 0.56537815, 0.52511588, 0.48429611,
               0.44291882, 0.50632691, 0.67452038, 0.        ] Jy>

    >>> spline = SplineInterpolatedResample()
    >>> new_spec3 = spline(input_spec, new_disp_grid)
    >>> new_spec3.flux #doctest:+SKIP
    <Quantity [0.18278657, 0.01050715, 0.27264392, 0.69624523, 1.01334485,
               1.07064075, 0.83599236, 0.43168126, 0.12365467, 0.05162207,
               0.21706343, 0.48931305, 0.62018392, 0.55961905, 0.45844228,
               0.4132248 , 0.45743191, 0.62178518, 0.        ] Jy>

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
    >>> spectral_axis = np.linspace(0., 10., 200) * u.GHz
    >>> spectral_model = models.Gaussian1D(amplitude=3*u.Jy, mean=5*u.GHz, stddev=0.8*u.GHz)
    >>> flux = spectral_model(spectral_axis)
    >>> flux += np.random.normal(0., 0.2, spectral_axis.shape) * u.Jy
    >>> noisy_gaussian = Spectrum1D(spectral_axis=spectral_axis, flux=flux)

Now we estimate the uncertainty from the region that does *not* contain
the line:

.. code-block:: python

    >>> from specutils import SpectralRegion
    >>> from specutils.manipulation import noise_region_uncertainty
    >>> noise_region = SpectralRegion([(0, 3), (7, 10)]*u.GHz)
    >>> spec_w_unc = noise_region_uncertainty(noisy_gaussian, noise_region)
    >>> spec_w_unc.uncertainty # doctest: +ELLIPSIS
    StdDevUncertainty([0.18461457, ..., 0.18461457])

Or similarly, expressed in pixels:

.. code-block:: python

    >>> noise_region = SpectralRegion([(0, 25), (175, 200)]*u.pix)
    >>> spec_w_unc = noise_region_uncertainty(noisy_gaussian, noise_region)
    >>> spec_w_unc.uncertainty # doctest: +ELLIPSIS
    StdDevUncertainty([0.18714535, ..., 0.18714535])

Reference/API
-------------

.. automodapi:: specutils.manipulation
    :no-heading:
