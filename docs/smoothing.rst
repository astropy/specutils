==================
Spectral Smoothing
==================

The spectral smoothing has two forms, 1) convolution based smoothing 
using the `astropy.convolution` package and 2) median filter
using the :func:`scipy.signal.medfilt`.  Each of these act on the flux
of the :class:`~specutils.spectra.Spectrum1D` object.

.. note:: Specutils smoothing kernel widths and standard deviations are
             in units of pixels and not ``Quantity``.

Convolution Based Smoothing
---------------------------

While any kernel supported by `astropy.convolution` will work (using the `~specutils.manipulation.convolution_smooth` function), 
several commonly-used kernels have convenience functions wrapping them to simplify the smoothing process into a simple 
one-line operation.  Currently implemented are:
:func:`~specutils.manipulation.box_smooth` (:class:`astropy.convolution.convolve.Box1DKernel`),  
:func:`~specutils.manipulation.gaussian_smooth` (:class:`astropy.convolution.convolve.Gaussian1DKernel`),  
and :func:`~specutils.manipulation.trapezoid_smooth` (:class:`astropy.convolution.convolve.Trapezoid1DKernel`).


.. code-block:: python

    >>> from specutils import Spectrum1D
    >>> import astropy.units as u
    >>> import numpy as np
    >>> from specutils.manipulation import (box_smooth, gaussian_smooth, trapezoid_smooth)

    >>> spec1 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm, flux=np.random.sample(49))
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

    >>> from specutils import Spectrum1D
    >>> import astropy.units as u
    >>> from astropy.convolution import Box1DKernel
    >>> import numpy as np
    >>> from specutils.manipulation import convolution_smooth

    >>> box1d_kernel = Box1DKernel(width=3)

    >>> spec1 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm, flux=np.random.sample(49))
    >>> spec1_bsmooth2 = convolution_smooth(spec1, box1d_kernel)

In this case, the ``spec1_bsmooth2`` result should be equivalent to the ``spec1_bsmooth`` in
the section above (assuming the flux data of the input ``spec`` is the same).


Median Smoothing
----------------

The median based smoothing  is implemented using `scipy.signal.medfilt` and
has a similar call structure to the convolution-based smoothing methods. This 
method applys the median filter across the flux.

Note: This method is not flux conserving.

.. code-block:: python

    >>> from specutils import Spectrum1D
    >>> import astropy.units as u
    >>> import numpy as np
    >>> from specutils.manipulation import median_smooth 

    >>> spec1 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm, flux=np.random.sample(49))
    >>> spec1_msmooth = median_smooth(spec1, width=3)

.. automodapi:: specutils.manipulation
