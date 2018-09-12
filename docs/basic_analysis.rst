.. currentmodule:: specutils.analysis

==============
Basic Analysis
==============

The specutils package comes with a few basic spectral analytic functions.
More extensive and involve analysis techniques will be available in another
package, `specreduce <https://github.com/astropy/specreduce>`_.

Specutils supports some built-in callable functions for basic calculations
over the given spectrum object.

Equivalent Width
----------------

Currently, specutils supports basic equivalent width calculations.

.. code-block:: python

    >>> import numpy as np
    >>> import astropy.units as u
    >>> from specutils.spectra import Spectrum1D
    >>> from specutils.analysis import equivalent_width

    >>> spec = Spectrum1D(spectral_axis=np.arange(50) * u.AA, flux=np.random.randn(50)*u.Jy)
    >>> equivalent_width(spec) #doctest:+SKIP
    <Quantity 24.16006697 Angstrom>

SNR
---

Currently, specutils supports basic signal-to-noise ratio calculations.

.. code-block:: python

    >>> import numpy as np
    >>> import astropy.units as u
    >>> from specutils.spectra import Spectrum1D
    >>> from astropy.nddata import StdDevUncertainty
    >>> from specutils.analysis import snr

    >>> uncertainty = StdDevUncertainty(0.1*np.abs(np.random.random(50))*u.Jy)
    >>> spec = Spectrum1D(spectral_axis=np.arange(50) * u.AA, flux=(3+np.random.randn(50))*u.Jy, uncertainty=uncertainty)
    >>> snr(spec) #doctest:+SKIP
    <Quantity 149.97247134>

Centroid
--------

Currently, specutils supports basic centroid calculations.

.. code-block:: python

    >>> import numpy as np
    >>> import astropy.units as u
    >>> from specutils.spectra import Spectrum1D
    >>> from astropy.nddata import StdDevUncertainty
    >>> from specutils.analysis import centroid

    >>> spec = Spectrum1D(spectral_axis=np.arange(50) * u.AA, flux=(3+np.random.randn(50))*u.Jy)
    >>> centroid(spec, region=None) #doctest:+SKIP
    <Quantity 24.39045495 Angstrom>

And if the spectrum contains a continuum, then it should be subtracted first:
.. code-block:: python

    >>> import numpy as np
    >>> import astropy.units as u
    >>> from specutils.spectra import Spectrum1D
    >>> from astropy.nddata import StdDevUncertainty
    >>> from specutils.fitting import fit_generic_continuum
    >>> from specutils.analysis import centroid

    >>> spec = Spectrum1D(spectral_axis=np.arange(50) * u.AA, flux=(10+np.random.randn(50))*u.Jy)
    >>> continuum_baseline = fit_generic_continuum(spec) #doctest:+SKIP
    >>> continuum_flux = continuum_baseline(spec.spectral_axis.value) #doctest:+SKIP
    >>> continuum = Spectrum1D(spectral_axis=spec.spectral_axis, flux=continuum_flux) #doctest:+SKIP
    >>> c = centroid(spec-continuum, region=None) #doctest:+SKIP

Width
-----

There are several width statistics that are provided by the
`~specutils.analysis` submodule.

The `~gaussian_sigma_width` function estimates the width of the spectrum by
computing an approximation of the standard deviation.

The `~gaussian_fwhm` function estimates the width of the spectrum at half max,
again by computing an approximation of the standard deviation.

Both of these functions assume that spectrum is approximately gaussian and that
it is centered on the spectral axis.

The function `~fwhm` provides an estimate of the full width of the spectrum at
half max that does not assume the spectrum is gaussian or centered on the
spectral axis. It locates the maximum, and then locates the value closest to
half of the maximum on either side, and measures the distance between them.

Consider the following noisy gaussian spectrum as an example:

.. plot::
   :include-source: true
   :context:

   >>> import numpy as np
   >>> from astropy import units as u
   >>> from astropy.modeling import models
   >>> from specutils import Spectrum1D
   >>> from specutils.analysis import gaussian_sigma_width
   >>> np.random.seed(0)

   >>> spectral_axis = np.linspace(0., 10., 200) * u.GHz
   >>> # Note that the gaussian is centered on the spectral axis
   >>> spectral_model = models.Gaussian1D(amplitude=3*u.Jy, mean=5*u.GHz, stddev=0.8*u.GHz)
   >>> flux = spectral_model(spectral_axis)
   >>> # Add noise
   >>> flux += np.random.normal(0., 0.2, spectral_axis.shape) * u.Jy
   >>> noisy_gaussian = Spectrum1D(spectral_axis=spectral_axis, flux=flux)

   >>> import matplotlib.pyplot as plt #doctest:+SKIP
   >>> plt.plot(noisy_gaussian.spectral_axis, noisy_gaussian.flux) #doctest:+SKIP

Each of the width analysis functions are applied to this spectrum below:

.. code-block:: python

   >>> from specutils.analysis import gaussian_sigma_width, gaussian_fwhm, fwhm
   >>> gaussian_sigma_width(noisy_gaussian)
   <Quantity 1.59661941 GHz>
   >>> gaussian_fwhm(noisy_gaussian)
   <Quantity 1.87987569 GHz>
   >>> fwhm(noisy_gaussian)
   <Quantity 1.85929648 GHz>

For uncentered spectra, `~fwhm` will generally provide a more accurate value
than either `~gaussian_sigma_width` or `~gaussian_fwhm`. Consider the following
example:

.. plot::
   :include-source: true
   :context:

   >>> np.random.seed(0)

   >>> spectral_axis = np.linspace(0., 10., 200) * u.GHz
   >>> # Note that the spectrum is not centered on the spectral axis
   >>> spectral_model = models.Gaussian1D(amplitude=3*u.Jy, mean=2*u.GHz, stddev=1.2*u.GHz)
   >>> flux = spectral_model(spectral_axis)
   >>> # Add noise
   >>> flux += np.random.normal(0., 0.2, spectral_axis.shape) * u.Jy
   >>> uncentered_gaussian = Spectrum1D(spectral_axis=spectral_axis, flux=flux)

   >>> plt.clf() #doctest:+SKIP
   >>> plt.plot(uncentered_gaussian.spectral_axis, uncentered_gaussian.flux) #doctest:+SKIP

Note that in this case, the `~fwhm` estimate is reasonable, whereas the others
are not:

.. code-block:: python

   >>> gaussian_sigma_width(uncentered_gaussian)
   <Quantity 5.23511731 GHz>
   >>> gaussian_fwhm(uncentered_gaussian)
   <Quantity 6.16387959 GHz>
   >>> fwhm(uncentered_gaussian)
   <Quantity 2.7638191 GHz>

Reference/API
-------------
.. automodapi:: specutils.analysis
    :no-heading:
