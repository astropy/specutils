================
Spectral Regions
================

A spectral region may be defined and may encompass one, or more,
sub-regions. They are defined independently of a `~specutils.Spectrum1D`
object in the sense that spectral regions like "near the Halpha line rest
wavelength" have meaning independent of the details of a particular spectrum.

Spectral regions can be defined either as a single region by passing
two `~astropy.units.Quantity`'s or by passing a list of 2-tuples. Note that
the units of these quantites can be any valid spectral unit *or* ``u.pixel``
(which indicates to use indexing directly).

.. code-block:: python

    >>> from astropy import units as u
    >>> from specutils.spectra import SpectralRegion

    >>> sr = SpectralRegion(0.45*u.um, 0.6*u.um)
    >>> sr_two = SpectralRegion([(0.45*u.um, 0.6*u.um), (0.8*u.um, 0.9*u.um)])

`~specutils.SpectralRegion` can be combined by using the '+' operator:

.. code-block:: python

    >>> from astropy import units as u
    >>> from specutils.spectra import SpectralRegion
    >>> sr = SpectralRegion(0.45*u.um, 0.6*u.um) + SpectralRegion(0.8*u.um, 0.9*u.um)

Regions can also be added in place:

.. code-block:: python

    >>> from astropy import units as u
    >>> from specutils.spectra import SpectralRegion

    >>> sr1 = SpectralRegion(0.45*u.um, 0.6*u.um)
    >>> sr2 = SpectralRegion(0.8*u.um, 0.9*u.um)
    >>> sr1 += sr2

Regions can be sliced by indexing by an integer or by a range:

.. code-block:: python

    >>> from astropy import units as u
    >>> from specutils.spectra import SpectralRegion

    >>> sr = SpectralRegion(0.15*u.um, 0.2*u.um) + SpectralRegion(0.3*u.um, 0.4*u.um) +\
    ...      SpectralRegion(0.45*u.um, 0.6*u.um) + SpectralRegion(0.8*u.um, 0.9*u.um) +\
    ...      SpectralRegion(1.0*u.um, 1.2*u.um) + SpectralRegion(1.3*u.um, 1.5*u.um)

    >>> # Get one spectral region (returns a SpectralRegion instance)
    >>> sone = sr1[0]

    >>> # Slice spectral region.
    >>> subsr = sr[3:5]
    >>> # SpectralRegion: 0.8 um - 0.9 um, 1.0 um - 1.2 um

The lower and upper bounds on a region are accessed by calling lower
or upper. The lower bound of a `~specutils.SpectralRegion` is the
minimum of the lower bounds of each sub-region and the upper bound is the
maximum of the upper bounds:

.. code-block:: python

    >>> from astropy import units as u
    >>> from specutils.spectra import SpectralRegion

    >>> sr = (SpectralRegion(0.15*u.um, 0.2*u.um) + SpectralRegion(0.3*u.um, 0.4*u.um) +
    ...       SpectralRegion(0.45*u.um, 0.6*u.um) + SpectralRegion(0.8*u.um, 0.9*u.um) +
    ...       SpectralRegion(1.0*u.um, 1.2*u.um) + SpectralRegion(1.3*u.um, 1.5*u.um))

    >>> # Bounds on the spectral region (most minimum and maximum bound)
    >>> sr.bounds
    (<Quantity 0.15 um>, <Quantity 1.5 um>)

    >>> # Lower bound on the spectral region (most minimum)
    >>> sr.lower
    <Quantity 0.15 um>

    >>> sr.upper
    <Quantity 1.5 um>

    >>> # Lower bound on one element of the spectral region.
    >>> sr[3].lower
    <Quantity 0.8 um>

One can also delete a sub-region:

.. code-block:: python

    >>> from astropy import units as u
    >>> from specutils.spectra import SpectralRegion

    >>> sr = (SpectralRegion(0.15*u.um, 0.2*u.um) + SpectralRegion(0.3*u.um, 0.4*u.um) +
    ...       SpectralRegion(0.45*u.um, 0.6*u.um) + SpectralRegion(0.8*u.um, 0.9*u.um) +
    ...       SpectralRegion(1.0*u.um, 1.2*u.um) + SpectralRegion(1.3*u.um, 1.5*u.um))

    >>> del sr[1]
    >>> sr
    Spectral Region, 5 sub-regions:
    (0.15 um, 0.2 um)   (0.45 um, 0.6 um)   (0.8 um, 0.9 um)
    (1.0 um, 1.2 um)    (1.3 um, 1.5 um)

There is also the ability to iterate:

.. code-block:: python

    >>> from astropy import units as u
    >>> from specutils.spectra import SpectralRegion

    >>> sr = (SpectralRegion(0.15*u.um, 0.2*u.um) + SpectralRegion(0.3*u.um, 0.4*u.um) +
    ...       SpectralRegion(0.45*u.um, 0.6*u.um) + SpectralRegion(0.8*u.um, 0.9*u.um) +
    ...       SpectralRegion(1.0*u.um, 1.2*u.um) + SpectralRegion(1.3*u.um, 1.5*u.um))

    >>> for s in sr:
    ...     print(s.lower)
    0.15 um
    0.3 um
    0.45 um
    0.8 um
    1.0 um
    1.3 um


And, lastly, there is the ability to invert a `~specutils.SpectralRegion` given a
lower and upper bound. For example, if a set of ranges are defined each defining a range
around lines, then calling invert will return a `~specutils.SpectralRegion` that
defines the baseline/noise regions:

.. code-block:: python

    >>> from astropy import units as u
    >>> from specutils.spectra import SpectralRegion

    >>> sr = (SpectralRegion(0.15*u.um, 0.2*u.um) + SpectralRegion(0.3*u.um, 0.4*u.um) +
    ...       SpectralRegion(0.45*u.um, 0.6*u.um) + SpectralRegion(0.8*u.um, 0.9*u.um) +
    ...       SpectralRegion(1.0*u.um, 1.2*u.um) + SpectralRegion(1.3*u.um, 1.5*u.um))

    >>> sr_inverted = sr.invert(0.05*u.um, 3*u.um)
    >>> sr_inverted
    Spectral Region, 7 sub-regions:
    (0.05 um, 0.15 um)   (0.2 um, 0.3 um)     (0.4 um, 0.45 um)
    (0.6 um, 0.8 um)     (0.9 um, 1.0 um)     (1.2 um, 1.3 um)
    (1.5 um, 3.0 um)

Region Extraction
-----------------

Given a `~specutils.SpectralRegion`, one can extract a sub-spectrum
from a `~specutils.Spectrum1D` object. If the `~specutils.SpectralRegion`
has multiple sub-regions then by default a list of `~specutils.Spectrum1D` objects
will be returned. If the ``return_single_spectrum`` argument is set to ``True``,
the resulting spectra will be concatenated together into a single
`~specutils.Spectrum1D` object instead.

An example of a single sub-region `~specutils.SpectralRegion`:

.. code-block:: python

    >>> from astropy import units as u
    >>> import numpy as np
    >>> from specutils import Spectrum1D, SpectralRegion
    >>> from specutils.manipulation import extract_region

    >>> region = SpectralRegion(8*u.nm, 22*u.nm)
    >>> spectrum = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm, flux=np.random.sample(49)*u.Jy)
    >>> sub_spectrum = extract_region(spectrum, region)
    >>> sub_spectrum.spectral_axis
    <SpectralAxis [ 8.,  9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21.,
       22.] nm>

Extraction also correctly interprets different kinds of spectral region units
as would be expected:

.. code-block:: python

    >>> from astropy import units as u
    >>> import numpy as np
    >>> from specutils import Spectrum1D, SpectralRegion
    >>> from specutils.manipulation import extract_region

    >>> spectrum = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm, flux=np.random.sample(49)*u.Jy)
    >>> region_angstroms = SpectralRegion(80*u.AA, 220*u.AA)
    >>> sub_spectrum = extract_region(spectrum, region_angstroms)
    >>> sub_spectrum.spectral_axis
    <SpectralAxis [ 8.,  9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21.,
       22.] nm>
    >>> region_pixels = SpectralRegion(7.5*u.pixel, 21.5*u.pixel)
    >>> sub_spectrum = extract_region(spectrum, region_pixels)
    >>> sub_spectrum.spectral_axis
    <SpectralAxis [ 8.,  9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21.,
       22.] nm>

An example of a multiple sub-region `~specutils.SpectralRegion`:

.. code-block:: python

    >>> from astropy import units as u
    >>> import numpy as np
    >>> from specutils import Spectrum1D, SpectralRegion
    >>> from specutils.manipulation import extract_region

    >>> region = SpectralRegion([(8*u.nm, 22*u.nm), (34*u.nm, 40*u.nm)])
    >>> spectrum = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm, flux=np.random.sample(49)*u.Jy)
    >>> sub_spectra = extract_region(spectrum, region)
    >>> sub_spectra[0].spectral_axis
    <SpectralAxis [ 8.,  9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21.,
       22.] nm>
    >>> sub_spectra[1].spectral_axis
    <SpectralAxis [34., 35., 36., 37., 38., 39., 40.] nm>

Multiple sub-regions can also be returned as a single concatenated spectrum:

.. code-block:: python

    >>> sub_spectrum = extract_region(spectrum, region, return_single_spectrum=True)
    >>> sub_spectrum.spectral_axis
    <SpectralAxis [ 8.,  9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21.,
        22., 34., 35., 36., 37., 38., 39., 40.] nm>

The bounding region that includes all data, including the ones that lie
in between disjointed spectral regions, can be extracted with
`specutils.manipulation.extract_bounding_spectral_region`:

.. code-block:: python

    >>> from astropy import units as u
    >>> import numpy as np
    >>> from specutils import Spectrum1D, SpectralRegion
    >>> from specutils.manipulation import extract_bounding_spectral_region

    >>> spectrum = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
    ...                       flux=np.random.default_rng(12345).random(49)*u.Jy)
    >>> region = SpectralRegion([(8*u.nm, 12*u.nm), (24*u.nm, 30*u.nm)])
    >>> sub_spectrum = extract_bounding_spectral_region(spectrum, region)
    >>> sub_spectrum.spectral_axis
    <SpectralAxis [ 8.,  9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21.,
        22., 23., 24., 25., 26., 27., 28., 29., 30.] nm>


`~specutils.manipulation.spectral_slab` is basically an alternate entry point for
`~specutils.manipulation.extract_region`. Notice the slightly different way to input the spectral
axis range to be extracted.
This function's purpose is to facilitate migration of ``spectral_cube`` functionality
into ``specutils``:

.. code-block:: python

    >>> from astropy import units as u
    >>> import numpy as np
    >>> from specutils import Spectrum1D, SpectralRegion
    >>> from specutils.manipulation import spectral_slab

    >>> spectrum = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
    ...                       flux=np.random.default_rng(12345).random(49)*u.Jy)
    >>> sub_spectrum = spectral_slab(spectrum, 8*u.nm, 20*u.nm)
    >>> sub_spectrum.spectral_axis
    <SpectralAxis [ 8.,  9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20.] nm>


Line List Compatibility
-----------------------

`~specutils.SpectralRegion` objects can also be created from the `~astropy.table.QTable` object returned from the line
finding functions:

.. code-block:: python

    >>> from astropy import units as u
    >>> import numpy as np
    >>> from specutils import Spectrum1D, SpectralRegion
    >>> from astropy.modeling.models import Gaussian1D
    >>> from specutils.fitting import find_lines_derivative

    >>> g1 = Gaussian1D(1, 4.6, 0.2)
    >>> g2 = Gaussian1D(2.5, 5.5, 0.1)
    >>> g3 = Gaussian1D(-1.7, 8.2, 0.1)

    >>> x = np.linspace(0, 10, 200)
    >>> y = g1(x) + g2(x) + g3(x)
    >>> spectrum = Spectrum1D(flux=y * u.Jy, spectral_axis=x * u.um)

    >>> lines = find_lines_derivative(spectrum, flux_threshold=0.01)
    >>> spec_reg = SpectralRegion.from_line_list(lines)
    >>> spec_reg
    Spectral Region, 3 sub-regions:
      (4.072864321608041 um, 5.072864321608041 um)
      (4.977386934673367 um, 5.977386934673367 um)
      (7.690954773869347 um, 8.690954773869347 um)

This can be fed into the ``exclude_regions`` argument of the `~specutils.fitting.fit_generic_continuum` or
`~specutils.fitting.fit_continuum` functions to avoid fitting regions that contain line features. Note that,
by default, this uses pythonic slicing, i.e., spectral values greater than or equal to the lower bound and
less than the upper bound of the region will be excluded from the fit. For convenience in some cases, the
``exclude_region_upper_bounds`` keyword can be set to ``True`` to exlude spectral values less than *or equal*
to the upper bound instead.

Conversely, users can also invert the spectral region

.. code-block:: python

    >>> inv_spec_reg = spec_reg.invert(spectrum.spectral_axis[0], spectrum.spectral_axis[-1])
    >>> inv_spec_reg
    Spectral Region, 3 sub-regions:
      (0.0 um, 4.072864321608041 um)
      (5.977386934673367 um, 7.690954773869347 um)
      (8.690954773869347 um, 10.0 um)

and use that result as the ``exclude_regions`` argument in the `~specutils.fitting.fit_lines` function in order to avoid
attempting to fit any of the continuum region.

Reading and Writing
-------------------

`~specutils.SpectralRegion` objects can be written to ``ecsv`` files, which uses the `~astropy.table.QTable` write machinery:

.. code-block:: python

    >>> spec_reg.write("spectral_region.ecsv")

This overwrites existing files by default if a duplicate filename is input. The resulting files can be read back in
to create a new `~specutils.SpectralRegion` using the ``read`` class method:

.. code-block:: python

    >>> spec_reg = SpectralRegion.read("spectral_region.ecsv")

.. testcleanup::

    >>> import os
    >>> os.remove("spectral_region.ecsv")

The `~astropy.table.QTable` created to write out the `~specutils.SpectralRegion` to file can also be accessed
directly with the ``as_table`` method, and a `~specutils.SpectralRegion` can be created directly from a `~astropy.table.QTable`
with the appropriate columns (minimally ``lower_bound`` and ``upper_bound``) using the ``from_qtable`` class method.

.. code-block:: python

    >>> spec_reg_table = spec_reg.as_table()
    >>> spec_reg_2 = SpectralRegion.from_qtable(spec_reg_table)

Reference/API
-------------

.. automodapi:: specutils
    :no-main-docstr:
    :no-heading:
    :no-inheritance-diagram:

    :skip: QTable
    :skip: Spectrum1D
    :skip: SpectrumCollection
    :skip: SpectralAxis
