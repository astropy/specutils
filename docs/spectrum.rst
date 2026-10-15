=============================
Working with Spectrum objects
=============================

As described in more detail in :doc:`types_of_spectra`, the core data class in
specutils for a single spectrum is :class:`~specutils.Spectrum`.  This object
can represent either one or many spectra, all with the same ``spectral_axis``.
This section describes some of the basic features of this class.

Basic Spectrum Creation
-----------------------

The simplest way to create a :class:`~specutils.Spectrum` is to
create it explicitly from arrays or :class:`~astropy.units.Quantity` objects:

.. plot::
    :include-source:
    :align: center

    >>> import numpy as np
    >>> import astropy.units as u
    >>> import matplotlib.pyplot as plt
    >>> from specutils import Spectrum
    >>> flux = np.random.randn(200)*u.Jy
    >>> wavelength = np.arange(5100, 5300)*u.AA
    >>> spec1d = Spectrum(spectral_axis=wavelength, flux=flux)
    >>> ax = plt.subplots()[1]  # doctest: +SKIP
    >>> ax.plot(spec1d.spectral_axis, spec1d.flux)  # doctest: +SKIP
    >>> ax.set_xlabel("Dispersion")  # doctest: +SKIP
    >>> ax.set_ylabel("Flux")  # doctest: +SKIP

.. note::
    The ``spectral_axis`` can also be provided as a :class:`~specutils.SpectralAxis` object,
    and in fact will internally convert the spectral_axis to :class:`~specutils.SpectralAxis` if it
    is provided as an array or `~astropy.units.Quantity`.

.. note::
    The ``spectral_axis`` can be either ascending or descending, but must be monotonic
    in either case.

Reading from a File
-------------------

``specutils`` takes advantage of the Astropy IO machinery and allows loading and
writing to files. The example below shows loading a FITS file.
``specutils`` has built-in (default) data loaders for some ASCII-based
formats and a range of FITS file formats specific to various
telescopes and observatories, but can also be extended by user's own
custom loaders (see below).

.. code-block:: python

    >>> from specutils import Spectrum
    >>> spec1d = Spectrum.read("/path/to/file.fits")  # doctest: +SKIP

Most of these default specutils loaders can also read an existing
`astropy.io.fits.HDUList` object (for FITS formats) or an open file object
(as resulting from e.g. streaming a file from the internet), and will
transparently support common compression formats such as ``gzip``,
``bzip2`` or ``lzma`` (``xz``).
Note that in these cases, a format string corresponding to an existing loader
should be supplied because these objects may lack enough contextual
information to automatically identify a loader.

.. code-block:: python

    >>> from specutils import Spectrum
    >>> import urllib
    >>> spec = urllib.request.urlopen('https://data.sdss.org/sas/dr14/sdss/spectro/redux/26/spectra/0751/spec-0751-52251-0160.fits') # doctest: +REMOTE_DATA
    >>> Spectrum.read(spec, format="SDSS-III/IV spec") # doctest: +REMOTE_DATA
    <Spectrum(flux=[30.59662628173828 ... 51.70271682739258] 1e-17 erg / (Angstrom s cm2) (shape=(3841,), mean=51.88042 1e-17 erg / (Angstrom s cm2)); spectral_axis=<SpectralAxis [3799.2686 3800.1426 3801.0188 ... 9193.905  9196.0205 9198.141 ] Angstrom> (length=3841); uncertainty=InverseVariance)>

Note that the same spectrum could be more conveniently downloaded via
astroquery, if the user has that package installed:

.. doctest-requires:: astroquery

     >>> from astroquery.sdss import SDSS  # doctest: +REMOTE_DATA
     >>> specs = SDSS.get_spectra(plate=751, mjd=52251, fiberID=160, data_release=14)  # doctest: +REMOTE_DATA
     >>> Spectrum.read(specs[0], format="SDSS-III/IV spec")  # doctest: +REMOTE_DATA
     <Spectrum(flux=[30.59662628173828 ... 51.70271682739258] 1e-17 erg / (Angstrom s cm2) (shape=(3841,), mean=51.88042 1e-17 erg / (Angstrom s cm2)); spectral_axis=<SpectralAxis [3799.2686 3800.1426 3801.0188 ... 9193.905  9196.0205 9198.141 ] Angstrom> (length=3841); uncertainty=InverseVariance)>


List of Loaders
~~~~~~~~~~~~~~~

The `~specutils.Spectrum` class has built-in support for various input and output formats.
A full list of the supported formats is shown in the table below and
can be accessed interactively with ``Spectrum.read.list_formats()``.
Note that the JWST readers require the ``stdatamodels`` package to be
installed, which is an optional dependency for ``specutils``.

.. automodule:: specutils.io._list_of_loaders

Call the help function for a specific loader to access further documentation
on that format and optional parameters accepted by the ``read`` function,
e.g. as ``Spectrum.read.help('tabular-fits')``. Additional optional parameters
are generally passed through to the backend functions performing the actual
reading operation, which depend on the loader. For loaders for FITS files for example,
this will often be :func:`astropy.io.fits.open`.

More information on creating custom loaders for formats not covered
by the above list can be found in the :doc:`custom loading </custom_loading>` page.

Writing to a File
-----------------

Similarly, a `~specutils.Spectrum` object can be saved to any of the
formats supporting writing (currently only the two generic FITS formats)
by using the :meth:`specutils.Spectrum.write` method.

.. code-block:: python

    >>> spec1d.write("/path/to/output.fits")  # doctest: +SKIP

Note that the above example, calling ``write()`` without specifying
any format, will default to the ``wcs1d-fits`` loader if the `~specutils.Spectrum`
has a compatible WCS, and to ``tabular-fits`` otherwise, or if writing
to another than the primary HDU (``hdu=0``) has been selected.
For better control of the file type, the ``format`` parameter should be explicitly passed.
Again, additional optional parameters are forwarded to the backend writing functions,
which for the FITS writers is :meth:`astropy.io.fits.HDUList.writeto`.

| More information on creating custom writers can be found in :ref:`custom_writer`.

Metadata
--------

The :attr:`specutils.Spectrum.meta` attribute provides a dictionary to store
additional information on the data, like origin, date and other circumstances.
For spectra read from files containing header-like attributes like a FITS
:class:`~astropy.io.fits.Header` or :attr:`astropy.table.Table.meta`,
loaders are conventionally storing this in ``Spectrum.meta['header']``.

The two provided FITS writers (``tabular-fits`` and ``wcs1d-fits``) save the contents of
``Spectrum.meta['header']`` (which should be an :class:`astropy.io.fits.Header`
or any object, like a `dict`, that can instantiate one) as the header of the
:class:`~astropy.io.fits.hdu.PrimaryHDU`.

Including Uncertainties
-----------------------

The :class:`~specutils.Spectrum` class supports
`Astropy uncertainties <https://docs.astropy.org/en/stable/nddata/nddata.html#uncertainties>`__.
Specifically, when initializing a :class:`~specutils.Spectrum`, any provided uncertainty should
be an instance of `~astropy.nddata.StdDevUncertainty`, `~astropy.nddata.VarianceUncertainty`,
or `~astropy.nddata.InverseVariance`. These subclasses of :class:`~astropy.nddata.NDData` have
propagation rules implemented at the class level, which are used by ``specutils`` to propagate
uncertainties when doing arithmetic on `~specutils.Spectrum` objects.

.. code-block:: python

    >>> from specutils import Spectrum
    >>> from astropy.nddata import StdDevUncertainty

    >>> spec = Spectrum(spectral_axis=np.arange(5000, 5010) * u.AA, flux=np.random.sample(10) * u.Jy, uncertainty=StdDevUncertainty(np.random.sample(10) * 0.1))

.. warning:: Not defining an uncertainty class will result in an
             :class:`~astropy.nddata.UnknownUncertainty` object which will not
             propagate uncertainties in arithmetic operations.


Including Masks
---------------

Masks are also available for :class:`~specutils.Spectrum`, following the
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
``mask`` attribute in :class:`~specutils.Spectrum` addresses that in that the
spectrum may still have a value underneath the mask, but it is not used in most
calculations. To allow for compatibility with ``NaN``-masking representations,
however, specutils will recognize ``flux`` values input as ``NaN`` and set the
mask to ``True`` for those values unless explicitly overridden.


Including Redshift or Radial Velocity
-------------------------------------

The :class:`~specutils.Spectrum` class supports setting a redshift or radial
velocity upon initialization of the object, as well as updating these values.
The default value for redshift and radial velocity is zero - to create a
:class:`~specutils.Spectrum` with a non-zero value, simply set the appropriate
attribute on object creation:

.. code-block:: python

    >>> spec1 = Spectrum(spectral_axis=np.arange(5000, 5010) * u.AA, flux=np.random.sample(10) * u.Jy, redshift = 0.15)
    >>> spec2 = Spectrum(spectral_axis=np.arange(5000, 5010) * u.AA, flux=np.random.sample(10) * u.Jy, radial_velocity = 1000 * u.Unit("km/s"))

By default, updating either the ``redshift`` or ``radial_velocity`` attributes
of an existing :class:`~specutils.Spectrum` directly uses the
:meth:`specutils.Spectrum.shift_spectrum_to` method, which also updates the
values of the ``spectral_axis`` to match the new frame. To leave the
``spectral_axis`` values unchanged while updating the ``redshift`` or
``radial_velocity`` value, use the :meth:`specutils.Spectrum.set_redshift_to`
or :meth:`specutils.Spectrum.set_radial_velocity_to` method as appropriate.
An example of the different treatments of the ``spectral_axis`` is shown below.

.. code-block:: python

    >>> spec1.shift_spectrum_to(redshift=0.5)  # Equivalent: spec1.redshift = 0.5
    >>> spec1.spectral_axis
    <SpectralAxis
       (observer to target:
          radial_velocity=115304.79153846155 km / s
          redshift=0.5000000000000002)
      [6521.73913043, 6523.04347826, 6524.34782609, 6525.65217391,
       6526.95652174, 6528.26086957, 6529.56521739, 6530.86956522,
       6532.17391304, 6533.47826087] Angstrom>
    >>> spec2.set_radial_velocity_to(5000 * u.Unit("km/s"))
    >>> spec2.spectral_axis
    <SpectralAxis
       (observer to target:
          radial_velocity=5000.0 km / s
          redshift=0.016819635148755285)
      [5000., 5001., 5002., 5003., 5004., 5005., 5006., 5007., 5008., 5009.] Angstrom>

Medium, Reference Frame and Observation Metadata
------------------------------------------------

The values of a ``spectral_axis`` only have a definite meaning once you know
the medium the wavelengths refer to (vacuum or air), the reference frame in
which they were measured (the telescope, the solar system barycentre, the
rest frame of the source, ...), and, to convert between frames, where the
source is and when and where it was observed. :class:`~specutils.Spectrum`
records all of these on its :class:`~specutils.SpectralAxis`, and they are
carried through slicing, arithmetic, resampling and region extraction:

.. code-block:: python

    >>> from astropy.coordinates import EarthLocation, SkyCoord
    >>> from astropy.time import Time
    >>> from specutils import SpectralMedium
    >>> apo = EarthLocation(lat=32.78 * u.deg, lon=-105.82 * u.deg, height=2788 * u.m)
    >>> spec = Spectrum(spectral_axis=np.linspace(5000, 5010, 11) * u.AA,
    ...                 flux=np.ones(11) * u.Jy,
    ...                 medium='air', frame='TOPOCENT',
    ...                 obstime=Time('2024-03-01T05:00:00'), location=apo,
    ...                 target=SkyCoord(ra=120 * u.deg, dec=-30 * u.deg))
    >>> spec.medium
    SpectralMedium(kind='air', refraction_method='Morton2000')
    >>> spec.frame
    'TOPOCENT'
    >>> spec.in_rest_frame
    False
    >>> spec[2:5].frame
    'TOPOCENT'

The metadata is:

* ``medium``: a :class:`~specutils.SpectralMedium`, or simply ``'vacuum'`` or
  ``'air'``. An air medium can also record which formula for the refractive
  index of air the wavelengths correspond to, and the air ``temperature``,
  ``pressure``, ``humidity`` and ``co2`` concentration (see
  `~specutils.utils.wcs_utils.refraction_index`). Air wavelengths cannot be
  converted to frequency, energy or velocity, since those relations only hold
  in vacuum; convert the spectrum to vacuum first.
* ``frame``: the reference frame the spectral values are measured in, using the
  FITS ``SPECSYS`` vocabulary: ``'TOPOCENT'``, ``'GEOCENTR'``, ``'BARYCENT'``,
  ``'HELIOCEN'``, ``'LSRK'``, ``'LSRD'``, ``'GALACTOC'``, ``'LOCALGRP'``,
  ``'CMBDIPOL'`` or ``'SOURCE'``, the rest frame of the source (see
  `~specutils.spectra.spectral_frame.SPECTRAL_FRAMES`; a few long forms such
  as ``'barycentric'`` or ``'rest'`` are also accepted). The ``radial_velocity``
  and ``redshift`` of the spectrum are the velocity of the source relative to an
  observer at rest in this frame, that is, the shift still to be applied to
  reach the ``'SOURCE'`` frame, where they are zero.
* ``target``: the position of the source as a
  :class:`~astropy.coordinates.SkyCoord`, optionally with its velocity.
* ``obstime`` and ``location``: the mid-point of the observation and the
  :class:`~astropy.coordinates.EarthLocation` of the telescope. Together with
  ``frame``, these define the ``observer`` whose rest frame the values are
  expressed in.

When both the observer and the target are known, the radial velocity is
computed from them, as for :class:`~astropy.coordinates.SpectralCoord`, so
``radial_velocity`` or ``redshift`` cannot also be given: set the velocity of
the source on the ``target`` instead. For a target with no velocity, the
topocentric radial velocity is just the motion of the observer, and the
``barycentric_correction`` is the velocity to add to velocities measured in
the current frame to obtain barycentric ones:

.. code-block:: python

    >>> spec.radial_velocity  # doctest: +FLOAT_CMP
    <Quantity 9.4102236 km / s>
    >>> spec.barycentric_correction  # doctest: +FLOAT_CMP
    <Quantity -9.4102236 km / s>

Converting between media and frames
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:meth:`~specutils.Spectrum.with_medium` returns a copy of the spectrum with
the spectral axis converted between vacuum and air (or between different air
conditions) using `~specutils.utils.wcs_utils.vac_to_air` and
`~specutils.utils.wcs_utils.air_to_vac`:

.. code-block:: python

    >>> vac = spec.with_medium('vacuum')
    >>> vac.medium
    SpectralMedium(kind='vacuum')
    >>> vac.spectral_axis.quantity[:3]  # doctest: +FLOAT_CMP
    <Quantity [5001.39486701, 5002.39513281, 5003.39539861] Angstrom>
    >>> thin_air = vac.with_medium(SpectralMedium('air', temperature=5 * u.deg_C,
    ...                                           pressure=700 * u.hPa))
    >>> thin_air.spectral_axis.quantity[:3]  # doctest: +FLOAT_CMP
    <Quantity [5000.39659203, 5001.3966676 , 5002.39674317] Angstrom>

:meth:`~specutils.Spectrum.with_frame` returns a copy transformed to another
reference frame, and :meth:`~specutils.Spectrum.to_rest` to the rest frame of
the source. When the observer and target are known the transformation is
computed from them; the result records its new ``frame`` and its
``radial_velocity`` is relative to that frame:

.. code-block:: python

    >>> bary = vac.with_frame('BARYCENT')
    >>> bary.frame, bary.radial_velocity, bary.barycentric_correction
    ('BARYCENT', <Quantity 0. km / s>, <Quantity 0. km / s>)
    >>> rest = vac.to_rest()
    >>> rest.frame, rest.in_rest_frame, rest.radial_velocity
    ('SOURCE', True, <Quantity 0. km / s>)

When they are not, the ``radial_velocity`` of the spectrum is enough to reach
the rest frame, and any other transformation needs the ``velocity`` to apply,
for instance the barycentric correction reported by a pipeline:

.. code-block:: python

    >>> obs = Spectrum(spectral_axis=np.linspace(5000, 5010, 11) * u.AA,
    ...                flux=np.ones(11) * u.Jy, frame='TOPOCENT',
    ...                radial_velocity=42 * u.km / u.s)
    >>> bary = obs.with_frame('BARYCENT', velocity=-9.41 * u.km / u.s)
    >>> bary.radial_velocity  # doctest: +FLOAT_CMP
    <Quantity 32.59 km / s>
    >>> bary.to_rest().radial_velocity
    <Quantity 0. km / s>

These methods replace the WCS of the spectrum with a lookup table built from
the new spectral axis, storing the original in ``meta['original_wcs']``, like
:meth:`~specutils.Spectrum.with_spectral_axis_unit`. Note that a velocity shift
of air wavelengths is only approximate (the refractive index is taken to be
constant over the shift), so convert to vacuum first, as above.

FITS keywords and loaders
~~~~~~~~~~~~~~~~~~~~~~~~~

A spectrum created from a FITS WCS takes its medium from the spectral
``CTYPE`` (``'AWAV'`` for air, ``'WAVE'``, ``'FREQ'``, ``'ENER'`` or
``'WAVN'`` for vacuum), its frame from ``SPECSYS``, its observation time from
``MJD-AVG`` (or the middle of ``MJD-BEG`` and ``MJD-END``, or ``MJD-OBS``) and
its location from the ``OBSGEO-*`` keywords. The ``wcs1d-fits`` and
``tabular-fits`` writers record these, along with the ``RA`` and ``DEC`` of
the target, and the corresponding loaders read them back.

Loaders for other formats can use
`~specutils.io.parsing_utils.spectral_axis_metadata_from_header`, which reads
the target position, the mid-point of the exposure and the location from the
usual header keywords, and takes what is known about the format (such as
``medium='vacuum'`` and ``frame='BARYCENT'``) explicitly:

.. code-block:: python

    >>> from specutils.io.parsing_utils import spectral_axis_metadata_from_header
    >>> header = {'RA': 120.0, 'DEC': -30.0, 'DATE-OBS': '2024-03-01T05:00:00',
    ...           'EXPTIME': 600.0}
    >>> metadata = spectral_axis_metadata_from_header(header, medium='vacuum',
    ...                                               frame='BARYCENT')
    >>> sorted(metadata)
    ['frame', 'medium', 'obstime', 'target']
    >>> metadata['obstime']
    <Time object: scale='utc' format='isot' value=2024-03-01T05:05:00.000>
    >>> spec = Spectrum(spectral_axis=np.linspace(5000, 5010, 11) * u.AA,
    ...                 flux=np.ones(11) * u.Jy, **metadata)
    >>> spec.frame
    'BARYCENT'

The metadata is also preserved when writing to and reading from ASDF.

.. _spectrum-defining-wcs:

Defining WCS
------------

Specutils always maintains a WCS object whether it is passed explicitly by the
user, or is created dynamically by specutils itself. In the latter case, the
user need not be aware that the WCS object is being used, and can interact
with the :class:`~specutils.Spectrum` object as if it were only a simple
data container.

Currently, specutils understands two WCS formats: FITS WCS and GWCS. When a user
does not explicitly supply a WCS object, specutils will fallback on an internal
GWCS object it will create.

.. note:: To create a custom adapter for a different WCS class (i.e. aside from
          FITSWCS or GWCS), please see the documentation on WCS Adapter classes.


Providing a FITS-style WCS
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    >>> from specutils.spectra import Spectrum
    >>> import astropy.wcs as fitswcs
    >>> import astropy.units as u
    >>> import numpy as np
    >>> my_wcs = fitswcs.WCS(header={
    ...     'CDELT1': 1, 'CRVAL1': 6562.8, 'CUNIT1': 'Angstrom', 'CTYPE1': 'WAVE',
    ...     'RESTFRQ': 1400000000, 'CRPIX1': 25})
    >>> spec = Spectrum(flux=[5,6,7] * u.Jy, wcs=my_wcs)
    >>> spec.spectral_axis  # doctest: +FLOAT_CMP
    <SpectralAxis
       (observer to target:
          radial_velocity=0.0 km / s
          redshift=0.0
        doppler_rest=1400000000.0 Hz
        doppler_convention=None)
      [6.5388e-07, 6.5398e-07, 6.5408e-07] m>
    >>> spec.wcs.pixel_to_world(np.arange(3))  # doctest: +FLOAT_CMP
    <SpectralCoord [6.5388e-07, 6.5398e-07, 6.5408e-07] m>

When creating a `~specutils.Spectrum` using a WCS, you can also use the
``move_spectral_axis`` argument to force the spectral axis to a certain dimension
of a multi-dimenasional flux array. Prior to ``specutils`` version 2.0, the flux
array was always reordered such that the spectral axis corresponded to the last
flux axis - this behavior can be reproduced by setting ``move_spectral_axis=-1``
or ``move_spectral_axis='last'``. Note that the relevant axes in the flux, mask,
and uncertainty arrays are simply swapped, and the swap is also reflected in the
resulting WCS. No check is currently done to ensure that the resulting array has
the spatial axes (most often RA and Dec) in any particular order.

Multi-dimensional Data Sets
---------------------------

`~specutils.Spectrum` also supports the multidimensional case where you
have, for example, an ``(n_spectra, n_pix)``
shaped data set where each ``n_spectra`` element provides a different flux
data array. ``flux`` and ``uncertainty`` may be multidimensional as
long as one dimension matches the shape of the spectral_axis. This is meant
to allow fast operations on collections of spectra that share the same
``spectral_axis``. While it may seem to conflict with the “1D” in the class
name, this name scheme is meant to communicate the presence of a single
common spectral axis. In cases where the flux axis corresponding to the spectral
axis cannot be determined automatically (for example, if multiple flux axes
have the same length as the spectral axis), the spectral axis must be specified
with the ``spectral_axis_index`` argument when initializing the
`~specutils.Spectrum`.

.. note:: The case where each flux data array is related to a *different* spectral
          axis is encapsulated in the :class:`~specutils.SpectrumCollection`
          object described in the :doc:`related docs </spectrum_collection>`.

.. code-block:: python

    >>> from specutils import Spectrum

    >>> spec = Spectrum(spectral_axis=np.arange(5000, 5010) * u.AA,
    ...                 flux=np.random.default_rng(12345).random((5, 10)) * u.Jy)
    >>> spec_slice = spec[0]
    >>> spec_slice.spectral_axis
    <SpectralAxis [5000., 5001., 5002., 5003., 5004., 5005., 5006., 5007., 5008., 5009.] Angstrom>
    >>> spec_slice.flux
    <Quantity [0.22733602, 0.31675834, 0.79736546, 0.67625467, 0.39110955,
               0.33281393, 0.59830875, 0.18673419, 0.67275604, 0.94180287] Jy>

While the above example only shows two dimensions, this concept generalizes to
any number of dimensions for `~specutils.Spectrum`.


Slicing
-------

As seen above, `~specutils.Spectrum` supports slicing in the same way as any
other array-like object. Additionally, a `~specutils.Spectrum` can be sliced
along the spectral axis using world coordinates.

.. code-block:: python

    >>> from specutils import Spectrum

    >>> spec = Spectrum(spectral_axis=np.arange(5000, 5010) * u.AA,
    ...                 flux=np.random.default_rng(12345).random((5, 10)) * u.Jy)
    >>> spec_slice = spec[5002*u.AA:5006*u.AA]
    >>> spec_slice.spectral_axis
    <SpectralAxis [5002., 5003., 5004., 5005.] Angstrom>

It is also possible to slice on other axes using simple array indices at the
same time as slicing the spectral axis based on spectral values.

.. code-block:: python

    >>> from specutils import Spectrum

    >>> spec = Spectrum(spectral_axis=np.arange(5000, 5010) * u.AA,
    ...                 flux=np.random.default_rng(12345).random((5, 10)) * u.Jy)
    >>> spec_slice = spec[2:4, 5002*u.AA:5006*u.AA]
    >>> spec_slice.shape
    (2, 4)

If the `specutils.Spectrum` was created with a WCS that included spatial
information, for example in case of a spectral cube with two spatial dimensions,
the `specutils.Spectrum.crop` method can be used to subset the data based on
the world coordinates. The inputs required are two sets up `astropy.coordinates`
objects defining the upper and lower corner of the region desired. Note that if
one of the coordinates is decreasing along an axis, the higher world coordinate
value will apply to the lower bound input.

.. code-block:: python

    >>> from astropy.coordinates import SpectralCoord, SkyCoord
    >>> from astropy import units as u
    >>> from astropy.wcs import WCS

    >>> w = WCS({'WCSAXES': 3, 'CRPIX1': 38.0, 'CRPIX2': 38.0, 'CRPIX3': 1.0,
    ...          'CRVAL1': 205.4384, 'CRVAL2': 27.004754, 'CRVAL3': 4.890499866509344,
    ...          'CTYPE1': 'RA---TAN', 'CTYPE2': 'DEC--TAN', 'CTYPE3': 'WAVE',
    ...          'CUNIT1': 'deg', 'CUNIT2': 'deg', 'CUNIT3': 'um',
    ...          'CDELT1': 3.61111097865634E-05, 'CDELT2': 3.61111097865634E-05, 'CDELT3': 0.001000000047497451,
    ...          'PC1_1 ': -1.0, 'PC1_2 ': 0.0, 'PC1_3 ': 0,
    ...          'PC2_1 ': 0.0, 'PC2_2 ': 1.0, 'PC2_3 ': 0,
    ...          'PC3_1 ': 0, 'PC3_2 ': 0, 'PC3_3 ': 1,
    ...          'DISPAXIS': 2, 'VELOSYS': -2538.02,
    ...          'SPECSYS': 'BARYCENT', 'RADESYS': 'ICRS', 'EQUINOX': 2000.0,
    ...          'LONPOLE': 180.0, 'LATPOLE': 27.004754})
    >>> spec = Spectrum(flux=np.random.default_rng(12345).random((20, 5, 10)) * u.Jy, wcs=w)  # doctest: +IGNORE_WARNINGS
    >>> lower = [SkyCoord(ra=205, dec=26, unit=u.deg), SpectralCoord(4.9, unit=u.um)]
    >>> upper = [SkyCoord(ra=205.5, dec=27.5, unit=u.deg), SpectralCoord(4.9, unit=u.um)]
    >>> spec.crop(lower, upper)  # doctest: +IGNORE_WARNINGS +FLOAT_CMP
    <Spectrum(flux=[[[0.708612359963129 ... 0.6345714580773677]]] Jy (shape=(1, 5, 10), mean=0.49653 Jy); spectral_axis=<SpectralAxis
        (observer to target:
           radial_velocity=0.0 km / s
           redshift=0.0)
      [4.90049987e-06] m> (length=1); medium=vacuum; frame=BARYCENT)>

Collapsing
----------

`~specutils.Spectrum` has built-in convenience methods for collapsing the
flux array of the spectrum via various statistics. The available statistics are
mean, median, sum, max, and min, and may be called either on a specific axis
(or axes) or over the entire flux array. The collapse methods currently respect
the ``mask`` attribute of the `~specutils.Spectrum`, but do not propagate
any ``uncertainty`` attached to the spectrum.

.. code-block:: python

    >>> spec = Spectrum(spectral_axis=np.arange(5000, 5010) * u.AA,
    ...                 flux=np.random.default_rng(12345).random((5, 10)) * u.Jy)
    >>> spec.mean()  # doctest: +FLOAT_CMP
    <Quantity 0.49802844 Jy>

The 'axis' argument of the collapse methods may either be an integer axis, or a
string specifying either 'spectral', which will collapse along only the
spectral axis, or 'spatial', which will collapse along all non-spectral axes.

.. code-block:: python

    >>> spec.mean(axis='spatial')  # doctest: +FLOAT_CMP
    <Spectrum(flux=<Quantity [0.37273938, 0.53843905, 0.61351648, 0.57311623, 0.44339915,
               0.66084728, 0.45881921, 0.38715911, 0.39967185, 0.53257671] Jy> (shape=(10,), mean=0.49803 Jy); spectral_axis=<SpectralAxis
      [5000. 5001. 5002. ... 5007. 5008. 5009.] Angstrom> (length=10))>

Note that in this case, the result of the collapse operation is a
`~specutils.Spectrum` rather than an `astropy.units.Quantity`, because the
collapse operation left the spectral axis intact.

It is also possible to supply your own function for the collapse operation by
calling `~specutils.Spectrum.collapse()` and providing a callable function
to the ``method`` argument.

.. code-block:: python

    >>> spec.collapse(method=np.nanmean, axis=1)  # doctest: +FLOAT_CMP
    <Quantity [0.51412398, 0.52665713, 0.31419772, 0.71556421, 0.41959918] Jy>

Reference/API
-------------

.. automodapi:: specutils
    :no-main-docstr:
    :inherited-members:
    :no-heading:
    :headings: -~

    :skip: SpectrumCollection
    :skip: SpectralRegion

.. automodapi:: specutils.spectra.spectral_frame
    :no-main-docstr:
    :no-heading:
    :no-inheritance-diagram:
    :include-all-objects:
    :headings: -~

    :skip: SpectralMedium
