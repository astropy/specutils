===================
Spectrum Collection
===================

A spectrum collection is a way to keep a set of :class:`~specutils.Spectrum1D`
objects together, regardless of their individual dispersion solutions. This is
done by resampling the spectra onto a defined output grid.

.. code:: python

    >>> from specutils import Spectrum1D, SpectrumCollection
    >>> import astropy.units as u
    >>> import numpy as np
    >>> spec1 = Spectrum1D(spectral_axis=np.linspace(0, 50, 50) * u.AA, flux=np.random.randn(50) * u.Jy)
    >>> spec2 = Spectrum1D(spectral_axis=np.linspace(0, 50, 25) * u.AA, flux=np.random.randn(25) * u.Jy)
    >>> spec_coll = SpectrumCollection([spec1, spec2], output_grid='coarse')
    >>> print(spec_coll.output_grid)

.. plot::

    import matplotlib.pyplot as plt
    from specutils import Spectrum1D, SpectrumCollection
    import astropy.units as u

    spec1 = Spectrum1D(spectral_axis=np.linspace(0, 50, 50) * u.AA,
                       flux=np.random.randn(50) * u.Jy)
    spec2 = Spectrum1D(spectral_axis=np.linspace(0, 50, 25) * u.AA,
                       flux=np.random.randn(25) * u.Jy)
    spec_coll = SpectrumCollection([spec1, spec2], output_grid='coarse')

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
    astropy.units.quantity.Quantity
    >>> print(type(spec_coll.flux))
    astropy.units.quantity.Quantity

The difference is their shape. The returned array from the :class:`~specutils.SpectrumCollection`
object will have shape ``(N, M)`` where ``N`` is the number of input spectra
and ``M`` is the length of the output dispersion grid.

.. code:: python

    >>> print(spec1.flux.shape)
    (25,)
    >>> print(spec_coll.flux.shape)
    (2, 24)

The items stored in the :class:`~specutils.SpectrumCollection` object are the
*original* input spectra. Iterating over a :class:`~specutils.SpectrumCollection`
will yield these original input spectra. Only when calling an attribute will
the user be returned a set of resampled values (or a set of values of length ``N``
in the case of e.g. ``WCS``). Therefore, to retrieve the resampled values of a
single spectrum, the user need only slice on the return array

.. code:: python

    >>> print(spec1.flux)  # Input spectrum values for flux
    >>> print(spec_coll[0].flux)  # Input spectrum values for flux
    >>> print(spec_coll.flux[0])  # Output spectrum values for flux


Defining output grids
---------------------

Output grids can be defined one of three ways: as a string indicating whether
to use the finest or coarsest sampling from the input spectra, or the to indicate
if all spectra should have the same dispersion grid (raising an exception if
this is not the case); an explicit array of values to use as the output; or as
a three-tuple indiciating the start bin, end bin, and bin size.

.. code:: python

    >>> spec_coll = SpectrumCollection([spec1, spec2], output_grid='fine')
    >>> print(spec_coll.wavelength.shape)
    (2, 50)
    >>> spec_coll = SpectrumCollection([spec1, spec2], output_grid='coarse')
    >>> print(spec_coll.wavelength.shape)
    (2, 24)
    >>> spec_coll = SpectrumCollection([spec1, spec2], output_grid=(0, 30, 1))
    >>> print(spec_coll.wavelength.shape)
    (2, 30)
    >>> spec_coll = SpectrumCollection([spec1, spec2], output_grid=[0., 1.42857143, 2.85714286, 4.28571429, 5.71428571, 7.14285714, 8.57142857, 10.])
    >>> print(spec_coll.wavelength.shape)
    (2, 8)

By default, the output grid is ``'same'``. This means that the
:class:`~specutils.SpectrumCollection` will maintain the bin size of the
spectra assuming they are all the same. Otherwise, an error will be raised
indicating that there is variation in the bin sizes.


Reference/API
-------------

.. automodapi:: specutils.spectra.spectrum_collection
    :no-heading:
