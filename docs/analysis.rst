.. currentmodule:: specutils.analysis

========
Analysis
========

The specutils package comes with a set of tools for doing common analysis
tasks on astronomical spectra. Some examples of applying these tools are
described below. The basic spectrum shown here is used in the examples in the
sub-sections below - a gaussian-profile line with flux of 5 GHz Jy.  See
:doc:`spectrum1d` for more on creating spectra:

.. plot::
    :include-source: true
    :context:

    >>> import numpy as np
    >>> from astropy import units as u
    >>> from astropy.nddata import StdDevUncertainty
    >>> from astropy.modeling import models
    >>> from specutils import Spectrum1D, SpectralRegion
    >>> np.random.seed(42)
    >>> spectral_axis = np.linspace(0., 10., 200) * u.GHz
    >>> spectral_model = models.Gaussian1D(amplitude=5*(2*np.pi*0.8**2)**-0.5*u.Jy, mean=5*u.GHz, stddev=0.8*u.GHz)
    >>> flux = spectral_model(spectral_axis)
    >>> flux += np.random.normal(0., 0.05, spectral_axis.shape) * u.Jy
    >>> uncertainty = StdDevUncertainty(0.2*np.ones(flux.shape)*u.Jy)
    >>> noisy_gaussian = Spectrum1D(spectral_axis=spectral_axis, flux=flux, uncertainty=uncertainty)
    >>> import matplotlib.pyplot as plt #doctest:+SKIP
    >>> plt.step(noisy_gaussian.spectral_axis, noisy_gaussian.flux) #doctest:+SKIP


SNR
---

The signal-to-noise ratio of a spectrum is often a valuable quantity for
evaluating the quality of a spectrum.  The `specutils.analysis.snr` function
performs this task, either on the spectrum as a whole, or sub-regions of a
spectrum:

.. code-block:: python

    >>> from specutils.analysis import snr
    >>> snr(noisy_gaussian)  # doctest:+FLOAT_CMP
    <Quantity 2.47730726>
    >>> snr(noisy_gaussian, SpectralRegion(4*u.GHz, 6*u.GHz))  # doctest:+FLOAT_CMP
    <Quantity 9.84136331>


A second method to calculate SNR does not require the uncertainty defined
on the `~specutils.Spectrum1D` object. This computes the signal to noise 
ratio DER_SNR following the definition set forth by the Spectral 
Container Working Group of ST-ECF, MAST and CADC. This is based on the
code at http://www.stecf.org/software/ASTROsoft/DER_SNR/.

.. code-block:: python

    >>> from specutils.analysis import snr_derived
    >>> snr_derived(noisy_gaussian)  # doctest:+FLOAT_CMP
    <Quantity 1.0448650809155666>
    >>> snr_derived(noisy_gaussian, SpectralRegion(4*u.GHz, 6*u.GHz))  # doctest:+FLOAT_CMP
    <Quantity 40.53528208761309>

The conditions on the data for this implementation for it to be an unbiased estimator of the SNR 
are strict.  In particular:

  * the noise is uncorrelated in wavelength bins spaced two pixels apart
  * for large wavelength regions, the signal over the scale of 5 or more pixels can be approximated by a straight line 


Line Flux Estimates
-------------------

While line-fitting (see :doc:`fitting`) is a more thorough way to measure
spectral line fluxes, direct measures of line flux are very useful for either
quick-look settings or for spectra not amedable to fitting.  The
`specutils.analysis.line_flux` function addresses that use case. The closely
related `specutils.analysis.equivalent_width` computes the equivalent width
of a spectral feature, a flux measure that is normalized against the continuum
of a spectrum.  Both are demonstrated below:

.. note::
    The `specutils.analysis.line_flux` function assumes the spectrum has
    already been continuum-subtracted, while
    `specutils.analysis.equivalent_width` assumes the continuum is at a fixed,
    known level (defaulting to 1, meaning continuum-normalized).
    :ref:`specutils-continuum-fitting` describes how continuua can be generated
    to prepare a spectrum for use with these functions.

.. code-block:: python

    >>> from specutils.analysis import line_flux
    >>> line_flux(noisy_gaussian).to(u.erg * u.cm**-2 * u.s**-1)  # doctest:+FLOAT_CMP
    <Quantity 4.97826405e-14 erg / (cm2 s)>
    >>> line_flux(noisy_gaussian, SpectralRegion(3*u.GHz, 7*u.GHz))  # doctest:+FLOAT_CMP
    <Quantity 4.92933252 GHz Jy>

For the equivalen width, note the need to add a continuum level:

.. code-block:: python

    >>> from specutils.analysis import equivalent_width
    >>> noisy_gaussian_with_continuum = noisy_gaussian + 1*u.Jy
    >>> equivalent_width(noisy_gaussian_with_continuum)  # doctest:+FLOAT_CMP
    <Quantity -4.97826405 GHz>
    >>> equivalent_width(noisy_gaussian_with_continuum, regions=SpectralRegion(3*u.GHz, 7*u.GHz))  # doctest:+FLOAT_CMP
    <Quantity -4.92933252 GHz>


Centroid
--------

The `specutils.analysis.centroid` function provides a first-moment analysis to
estimate the center of a spectral feature:

.. code-block:: python

    >>> from specutils.analysis import centroid
    >>> centroid(noisy_gaussian, SpectralRegion(3*u.GHz, 7*u.GHz))  # doctest:+FLOAT_CMP
    <Quantity 4.99881315 GHz>

While this example is "pre-subtracted", this function only performs well if the
contiuum has already been subtracted, as for the other functions above and
below.

Line Widths
-----------

There are several width statistics that are provided by the
`specutils.analysis` submodule.

The `~gaussian_sigma_width` function estimates the width of the spectrum by
computing a second-moment-based approximation of the standard deviation.

The `~gaussian_fwhm` function estimates the width of the spectrum at half max,
again by computing an approximation of the standard deviation.

Both of these functions assume that the spectrum is approximately gaussian.

The function `~fwhm` provides an estimate of the full width of the spectrum at
half max that does not assume the spectrum is gaussian. It locates the maximum,
and then locates the value closest to half of the maximum on either side, and
measures the distance between them.

A function to calculate the full width at zero intensity (i.e. the width of a
spectral feature at the continuum) is provided as `~fwzi`. Like the `~fwhm`
calculation, it does not make assumptions about the shape of the feature
and calculates the width by finding the points at either side of maximum
that reach the continuum value. In this case, it assumes the provided
spectrum has been continuum subtracted.

Each of the width analysis functions are applied to this spectrum below:

.. code-block:: python

   >>> from specutils.analysis import gaussian_sigma_width, gaussian_fwhm, fwhm, fwzi
   >>> gaussian_sigma_width(noisy_gaussian) # doctest: +FLOAT_CMP
   <Quantity 0.76925064 GHz>
   >>> gaussian_fwhm(noisy_gaussian) # doctest: +FLOAT_CMP
   <Quantity 1.81144683 GHz>
   >>> fwhm(noisy_gaussian) # doctest: +FLOAT_CMP
   <Quantity 1.91686888 GHz>
   >>> fwzi(noisy_gaussian) # doctest: +FLOAT_CMP
   <Quantity 89.99999455 GHz>


Template Comparison
-------------------

With one observed spectrum and n template spectra, the process to do template matching is:

    1. Move the templates as close as possible to the observed spectrum by doing redshifting, matching the resolution,
       and matching the wavelength spacing
    2. Then loop through all of the n template spectra and run chi square on them using the observed spectrum as the
       expected frequency and the template spectrum as the observed
    3. Once you have a corresponding chi square for each template spectrum, return the lowest chi square and its
       corresponding template spectrum (normalized) and the index of the template spectrum if the original template
       parameter is iterable

An example of how to do template matching in is:

.. code-block:: python

   >>> from specutils.analysis import template_comparison
   >>> spec_axis = np.linspace(0, 50, 50) * u.AA
   >>> observed_spectrum = Spectrum1D(spectral_axis=spec_axis, flux=np.random.randn(50) * u.Jy, uncertainty=StdDevUncertainty(np.random.sample(50), unit='Jy'))
   >>> spectral_template = Spectrum1D(spectral_axis=spec_axis, flux=np.random.randn(50) * u.Jy, uncertainty=StdDevUncertainty(np.random.sample(50), unit='Jy'))
   >>> tm_result = template_comparison.template_match(observed_spectrum, spectral_template) # doctest:+FLOAT_CMP


Reference/API
-------------
.. automodapi:: specutils.analysis
    :no-heading:
