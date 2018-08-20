===================
Spectrum Collection
===================

A spectrum collection is a way to keep a set of :class:`~specutils.Spectrum1D`
objects together and have the collection behave as if it were a single spectrum
object. This means that it can be used in regular analysis functions to perform
operations over entire sets of data.

Currently, all :class:`~specutils.SpectrumCollection` items must be the same
shape. No assumptions are made about the dispersion solutions, and users are
encouraged to ensure their spectrum collections make sense either by resampling
them beforehand, or being aware that they do not defautly share the same
dispersion solution.

.. code:: python

    >>> from specutils import Spectrum1D, SpectrumCollection
    >>> import astropy.units as u
    >>> import numpy as np
    >>> spec1 = Spectrum1D(spectral_axis=np.linspace(0, 50, 50) * u.AA, flux=np.random.randn(50) * u.Jy)
    >>> spec2 = Spectrum1D(spectral_axis=np.linspace(25, 76, 50) * u.AA, flux=np.random.randn(50) * u.Jy)
    >>> spec_coll = SpectrumCollection([spec1, spec2])
    >>> print(spec_coll.shape)
    [(50, 50), (50, 50)]

.. plot::

    import matplotlib.pyplot as plt
    from specutils import Spectrum1D, SpectrumCollection
    import astropy.units as u

    spec1 = Spectrum1D(spectral_axis=np.linspace(0, 50, 50) * u.AA,
                       flux=np.random.randn(50) * u.Jy)
    spec2 = Spectrum1D(spectral_axis=np.linspace(25, 76, 50) * u.AA,
                       flux=np.random.randn(50) * u.Jy)
    spec_coll = SpectrumCollection([spec1, spec2])

    f, (ax1, ax2) = plt.subplots(2, 1)

    ax1.plot(spec1.wavelength, spec1.flux, label="High-res")
    ax1.plot(spec2.wavelength, spec2.flux, label="Low-res")
    ax2.plot(spec_coll.wavelength[0], spec_coll.flux[0])
    ax2.plot(spec_coll.wavelength[1], spec_coll.flux[1])

:class:`~specutils.SpectrumCollection` objects can be treated just like
:class:`~specutils.Spectrum1D` objects; calling a particular attribute on the
object will return an array whose type depends on the type of the attribute in
the :class:`~specutils.Spectrum1D` object.

.. code:: python

    >>> print(type(spec1.flux))
    <class 'astropy.units.quantity.Quantity'>
    >>> print(type(spec_coll.flux))
    <class 'astropy.units.quantity.Quantity'>

The difference is their shape. The returned array from the
:class:`~specutils.SpectrumCollection`
object will have shape ``(N, M)`` where ``N`` is the number of input spectra
and ``M`` is the length of the output dispersion grid.

.. code:: python

    >>> print(spec1.flux.shape)
    (50,)
    >>> print(spec_coll.flux.shape)
    (2, 50)

The items stored in the :class:`~specutils.SpectrumCollection` object are the
*original* input spectra. Iterating over a :class:`~specutils.SpectrumCollection`
will yield these original input spectra. Only when calling an attribute will
the user be returned a set of values -- either numpy arrays, or lists of objects
(e.g. in the case of the spectra's wcs).


Reference/API
-------------

.. automodapi:: specutils.spectra.spectrum_collection
    :no-heading:
