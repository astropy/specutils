==================
Spectral Smoothing
==================

The spectral smoothing has two forms, 1) convolution based smoothing 
using the ``astropy.convolution`` package and 2) median filter
using the :func:``scipy.signal.medfilt``.  Each of these act on the flux
of the :class:`~specutils.spectra.Spectram1D` object.

Convolution Based Smoothing
---------------------------

The convolution based smoothing supports smoothing based on 1D kernels in 
``astropy.convlution.convolve``.  Currently implemented are 
:func:`~specutils.processing.smoothing.box_smooth` (:class:``astropy.convolution.convolve.Box1DKernel``),  
:func:`~specutils.processing.smoothing.gaussian_smooth` (:class:``astropy.convolution.convolve.Gaussian1DKernel``),  
and :func:`~specutils.processing.smoothing.trapzoid_smooth` (:class:``astropy.convolution.convolve.Trapezoid1DKernel``).

Each of the functions also has a parameter ``inplace`` which takes a boolean 
to tell the method to either apply the filter to the flux on the input object
(``inplace=True``) or to filter the flux on a copy of the input object (``inplace=False``).
In the latter case, the returned spectrum is a copy of the input spectrum but with the
smoothed flux.

.. code-block:: python

    >>> from specutils import Spectrum1D
    >>> import astropy.units as u
    >>> import numpy as np
    >>> from specutils.processing.smoothing import (box_smooth, 
                               gaussian_smooth, trapezoid_smooth)

    >>> spec1 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm, flux=np.random.sample(49))
    >>> spec1_bsmooth = box_smooth(spec1, width=3)
    >>> spec1_gsmooth = gaussian_smooth(spec1, stddev=3)
    >>> spec1_tsmooth = trapezoid_smooth(spec1, width=3)

    >>> gaussian_smooth(spec1, width=3, inplace=True)

Each of the specific smoothing methods create the appropriate ``astropy.convolution.convolve`` 
kernel and then call a helper function :func:`~specutils.processing.smoothing.convolution_smooth` 
that takes the spectrum and an astropy 1D kernel.  So, one could also do:

.. code-block:: python

    >>> from specutils import Spectrum1D
    >>> import astropy.units as u
    >>> from astropy.convolution.convolve import Box1DKernel
    >>> import numpy as np
    >>> from specutils.processing.smoothing import convolution_smooth

    >>> box1d_kernel = convolution.Box1DKernel(width)

    >>> spec1 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm, flux=np.random.sample(49))
    >>> spec1_bsmooth2 = convolution_smooth(spectrum, box1d_kernel, inplace)

In this case, the ``spec1_bsmooth2`` result should be equivalent to the ``spec1_bsmooth`` in
the section above (assuming the flux data of the input `spec` is the same).


Median Smoothing
----------------

The median based smoothing  is implemented using ``scipy.signal.medfilt`` and
has a similar call structure to the convolution-based smoothing methods. This 
method applys the median filter across the flux.

Note: This method is not flux conserving.

.. code-block:: python

    >>> from specutils import Spectrum1D
    >>> import astropy.units as u
    >>> import numpy as np
    >>> from specutils.processing.smoothing import median_smooth 

    >>> spec1 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm, flux=np.random.sample(49))
    >>> spec1_msmooth = median_smooth(spec1, width=3)
