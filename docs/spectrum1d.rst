========================
Working with Spectrum1Ds
========================

As described in more detail in :doc:`types_of_spectra`, the core data class in
specutils for a single spectrum is `~specutils.Spectrum1D`.  This object
can represent either one or many spectra, all with the same ``spectral_axis``.
This section describes some of the basic features of this class.

Basic Spectrum Creation
-----------------------

The simplest (and most powerful) way to create a `~specutils.Spectrum1D` is to
create it explicitly from arrays or `~astropy.units.Quantity` objects:

.. plot::
    :include-source:
    :align: center

    >>> import numpy as np
    >>> import astropy.units as u
    >>> import matplotlib.pyplot as plt
    >>> from specutils import Spectrum1D
    >>> flux = np.random.randn(200)*u.Jy
    >>> wavelength = np.arange(5100, 5300)*u.AA
    >>> spec1d = Spectrum1D(spectral_axis=wavelength, flux=flux)
    >>> ax = plt.subplots()[1]  # doctest: +SKIP
    >>> ax.plot(spec1d.spectral_axis, spec1d.flux)  # doctest: +SKIP
    >>> ax.set_xlabel("Dispersion")  # doctest: +SKIP
    >>> ax.set_ylabel("Flux")  # doctest: +SKIP

.. note::
    Note that the ``spectral_axis`` can also be provided as a :class:`~specutils.SpectralAxis` object,
    and in fact will internally convert the spectral_axis to :class:`~specutils.SpectralAxis` if it
    is provided as an array or `~astropy.units.Quantity`.

Reading from a File
-------------------

``specutils`` takes advantage of the Astropy IO machinery and allows loading and
writing to files. The example below shows loading a FITS file. While specutils
has some basic data loaders, for more complicated or custom files, users are
encouraged to :doc:`create their own loader </custom_loading>`.

.. code-block:: python

    >>> from specutils import Spectrum1D
    >>> spec1d = Spectrum1D.read("/path/to/file.fits")  # doctest: +SKIP

Most of the built-in specutils default loaders can also read an existing
`astropy.io.fits.HDUList` object or an open file object (as resulting
from e.g. streaming a file from the internet). Note that in these cases, a
format string corresponding to an existing loader must be supplied because
these objects lack enough contextual information to automatically identify
a loader.

.. code-block:: python

    >>> from specutils import Spectrum1D
    >>> import urllib
    >>> specs = urllib.request.urlopen('https://data.sdss.org/sas/dr14/sdss/spectro/redux/26/spectra/0751/spec-0751-52251-0160.fits') # doctest: +REMOTE_DATA
    >>> Spectrum1D.read(specs, format="SDSS-III/IV spec") # doctest: +REMOTE_DATA
    <Spectrum1D(flux=<Quantity [30.596626,...]...>

Note that the same spectrum could be more conveniently downloaded via
astroquery, if the user has that package installed:

.. code-block:: python

     >>> from astroquery.sdss import SDSS # doctest: +SKIP
     >>> specs = SDSS.get_spectra(plate=751, mjd=52251, fiberID=160) # doctest: +SKIP
     >>> Spectrum1D.read(specs[0], format="SDSS-III/IV spec") # doctest: +SKIP
     <Spectrum1D(flux=<Quantity [30.596626,...]...>


List of Loaders
~~~~~~~~~~~~~~~

The `~specutils.Spectrum1D` class has built-in support for various input and output formats.
A full list of the supported formats is shown in the table below.

.. automodule:: specutils.io._list_of_loaders

| More information on creating custom loaders can be found in the :doc:`custom loading </custom_loading>` page.


Including Uncertainties
-----------------------

The :class:`~specutils.Spectrum1D` class supports uncertainties, and
arithmetic operations performed with :class:`~specutils.Spectrum1D`
objects will propagate uncertainties.

Uncertainties are a special subclass of :class:`~astropy.nddata.NDData`, and their
propagation rules are implemented at the class level. Therefore, users must
specify the uncertainty type at creation time

.. code-block:: python

    >>> from specutils import Spectrum1D
    >>> from astropy.nddata import StdDevUncertainty

    >>> spec = Spectrum1D(spectral_axis=np.arange(5000, 5010)*u.AA, flux=np.random.sample(10)*u.Jy, uncertainty=StdDevUncertainty(np.random.sample(10) * 0.1))

.. warning:: Not defining an uncertainty class will result in an
             :class:`~astropy.nddata.UnknownUncertainty` object which will not
             propagate uncertainties in arithmetic operations.


Including Masks
---------------

Masks are also available for :class:`~specutils.Spectrum1D`, following the
same mechanisms as :class:`~astropy.nddata.NDData`.  That is, the mask should
have the property that it is ``False``/``0`` wherever the data is *good*, and
``True``/anything else where it should be masked.  This allows "data quality"
arrays to function as masks by default.

Note that this is distinct from "missing data" implementations, which generally
use ``NaN`` as a masking technique.  This method has the problem that ``NaN``
values are frequently "infectious", in that arithmetic operations sometimes
propagate to yield results as just ``NaN`` where the intent is instead to skip
that particular pixel. It also makes it impossible to store data that in the
spectrum that may have meaning but should *sometimes* be masked.  The separate
``mask`` attribute in :class:`~specutils.Spectrum1D` addresses that in that the
spectrum may still have a value underneath the mask, but it is not used in most
calculations. To allow for compatibility with ``NaN``-masking representations,
however, specutils will recognize ``flux`` values input as ``NaN`` and set the
mask to ``True`` for those values unless explicitly overridden.


Defining WCS
------------

Specutils always maintains a WCS object whether it is passed explicitly by the
user, or is created dynamically by specutils itself. In the latter case, the
user need not be aware that the WCS object is being used, and can interact
with the :class:`~specutils.Spectrum1D` object as if it were only a simple
data container.

Currently, specutils understands two WCS formats: FITS WCS and GWCS. When a user
does not explicitly supply a WCS object, specutils will fallback on an internal
GWCS object it will create.

.. note:: To create a custom adapter for a different WCS class (i.e. aside from
          FITSWCS or GWCS), please see the documentation on WCS Adapter classes.


Providing a FITS-style WCS
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    >>> from specutils.spectra import Spectrum1D
    >>> import astropy.wcs as fitswcs
    >>> import astropy.units as u
    >>> import numpy as np
    >>> my_wcs = fitswcs.WCS(header={'CDELT1': 1, 'CRVAL1': 6562.8, 'CUNIT1': 'Angstrom', 'CTYPE1': 'WAVE', 'RESTFRQ': 1400000000, 'CRPIX1': 25})
    >>> spec = Spectrum1D(flux=[5,6,7] * u.Jy, wcs=my_wcs)
    >>> spec.spectral_axis #doctest:+SKIP
    <Quantity [ 6538.8, 6539.8, 6540.8] Angstrom>
    >>> spec.wcs.pixel_to_world(np.arange(3)) #doctest:+SKIP
    array([6.5388e-07, 6.5398e-07, 6.5408e-07])


Multi-dimensional Data Sets
---------------------------

`~specutils.Spectrum1D` also supports the multidimensional case where you
have, say, an ``(n_spectra, n_pix)``
shaped data set where each ``n_spectra`` element provides a different flux
data array and so ``flux`` and ``uncertainty`` may be multidimensional as
long as the last dimension matches the shape of spectral_axis This is meant
to allow fast operations on collections of spectra that share the same
``spectral_axis``. While it may seem to conflict with the “1D” in the class
name, this name scheme is meant to communicate the presence of a single
common spectral axis.

.. note:: The case where each flux data array is related to a *different* spectral
          axis is encapsulated in the :class:`~specutils.SpectrumCollection`
          object described in the :doc:`related docs </spectrum_collection>`.

.. code-block:: python

    >>> from specutils import Spectrum1D

    >>> spec = Spectrum1D(spectral_axis=np.arange(5000, 5010)*u.AA, flux=np.random.sample((5, 10))*u.Jy)
    >>> spec_slice = spec[0] #doctest:+SKIP
    >>> spec_slice.spectral_axis #doctest:+SKIP
    <Quantity [5000., 5001., 5002., 5003., 5004., 5005., 5006., 5007., 5008., 5009.] Angstrom>
    >>> spec_slice.flux #doctest:+SKIP
    <Quantity [0.72722821, 0.32147784, 0.70256482, 0.04445197, 0.03390352,
           0.50835299, 0.87581725, 0.50270413, 0.08556376, 0.53713355] Jy>

While the above example only shows two dimensions, this concept generalizes to
any number of dimensions for `~specutils.Spectrum1D`, as long as the spectral
axis is always the last.



Reference/API
-------------

.. automodapi:: specutils
    :no-main-docstr:
    :inherited-members:
    :no-heading:
    :headings: -~

    :skip: test
    :skip: SpectrumCollection
    :skip: SpectralRegion
    :skip: UnsupportedPythonError
