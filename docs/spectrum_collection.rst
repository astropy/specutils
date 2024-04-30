================================
Working With SpectrumCollections
================================

A spectrum collection is a way to keep a set of spectra data together and have
the collection behave as if it were a single spectrum object. This means that
it can be used in regular analysis functions to perform operations over entire
sets of data.

Currently, all :class:`~specutils.SpectrumCollection` items must be the same
shape. No assumptions are made about the dispersion solutions, and users are
encouraged to ensure their spectrum collections make sense either by resampling
them beforehand, or being aware that they do not share the same dispersion
solution.

.. code:: python

    >>> import numpy as np
    >>> import astropy.units as u
    >>> from astropy.nddata import StdDevUncertainty
    >>> from specutils import SpectrumCollection
    >>> from specutils.utils.wcs_utils import gwcs_from_array

    >>> flux = u.Quantity(np.random.sample((5, 10)), unit='Jy')
    >>> spectral_axis = u.Quantity(np.arange(50).reshape((5, 10)), unit='AA')
    >>> wcs = np.array([gwcs_from_array(x) for x in spectral_axis])
    >>> uncertainty = StdDevUncertainty(np.random.sample((5, 10)), unit='Jy')
    >>> mask = np.ones((5, 10)).astype(bool)
    >>> meta = [{'test': 5, 'info': [1, 2, 3]} for i in range(5)]

    >>> spec_coll = SpectrumCollection(
    ... flux=flux, spectral_axis=spectral_axis, wcs=wcs,
    ... uncertainty=uncertainty, mask=mask, meta=meta)

    >>> spec_coll.shape
    (5,)
    >>> spec_coll.flux.unit
    Unit("Jy")
    >>> spec_coll.spectral_axis.shape
    (5, 10)
    >>> spec_coll.spectral_axis.unit
    Unit("Angstrom")

Collections from 1D spectra
---------------------------

It is also possible to create a :class:`~specutils.SpectrumCollection` from
a list of :class:`~specutils.Spectrum1D`:

.. code:: python

    >>> import warnings
    >>> import numpy as np
    >>> from astropy import units as u
    >>> from specutils import Spectrum1D, SpectrumCollection
    >>> spec = Spectrum1D(spectral_axis=np.linspace(0, 50, 50) * u.AA,
    ...                   flux=np.random.randn(50) * u.Jy,
    ...                   uncertainty=StdDevUncertainty(np.random.sample(50), unit='Jy'))
    >>> with warnings.catch_warnings():  # Ignore warnings
    ...     warnings.simplefilter('ignore')
    ...     spec1 = Spectrum1D(spectral_axis=np.linspace(20, 60, 50) * u.AA,
    ...                        flux=np.random.randn(50) * u.Jy,
    ...                        uncertainty=StdDevUncertainty(np.random.sample(50), unit='Jy'))
    ...     spec_coll = SpectrumCollection.from_spectra([spec, spec1])

    >>> spec_coll.shape
    (2,)
    >>> spec_coll.flux.unit
    Unit("Jy")
    >>> spec_coll.spectral_axis.shape
    (2, 50)
    >>> spec_coll.spectral_axis.unit
    Unit("Angstrom")

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


Reference/API
-------------

.. automodapi:: specutils
    :no-main-docstr:
    :no-heading:
    :no-inheritance-diagram:

    :skip: test
    :skip: QTable
    :skip: Spectrum1D
    :skip: SpectralRegion
    :skip: SpectralAxis
