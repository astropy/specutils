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

    >>> # Get on spectral region (returns a SpectralRegion instance)
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

    >>> sr = SpectralRegion(0.15*u.um, 0.2*u.um) + SpectralRegion(0.3*u.um, 0.4*u.um) +\
    ...      SpectralRegion(0.45*u.um, 0.6*u.um) + SpectralRegion(0.8*u.um, 0.9*u.um) +\
    ...      SpectralRegion(1.0*u.um, 1.2*u.um) + SpectralRegion(1.3*u.um, 1.5*u.um)

    >>> # Bounds on the spectral region (most minimum and maximum bound)
    >>> print(sr.bounds) #doctest:+SKIP
    (<Quantity 0.15 um>, <Quantity 1.5 um>)

    >>> # Lower bound on the spectral region (most minimum)
    >>> sr.lower #doctest:+SKIP
    <Quantity 0.15 um>

    >>> sr.upper #doctest:+SKIP
    <Quantity 1.5 um>

    >>> # Lower bound on one element of the spectral region.
    >>> sr[3].lower #doctest:+SKIP
    <Quantity 0.8 um>

One can also delete a sub-region:

.. code-block:: python

    >>> from astropy import units as u
    >>> from specutils.spectra import SpectralRegion

    >>> sr = SpectralRegion(0.15*u.um, 0.2*u.um) + SpectralRegion(0.3*u.um, 0.4*u.um) +\
    ...      SpectralRegion(0.45*u.um, 0.6*u.um) + SpectralRegion(0.8*u.um, 0.9*u.um) +\
    ...      SpectralRegion(1.0*u.um, 1.2*u.um) + SpectralRegion(1.3*u.um, 1.5*u.um)

    >>> del sr[1]
    >>> sr #doctest:+SKIP
    Spectral Region, 5 sub-regions:
    (0.15 um, 0.2 um)   (0.45 um, 0.6 um)   (0.8 um, 0.9 um)
    (1.0 um, 1.2 um)    (1.3 um, 1.5 um)

There is also the ability to iterate:

.. code-block:: python

    >>> from astropy import units as u
    >>> from specutils.spectra import SpectralRegion

    >>> sr = SpectralRegion(0.15*u.um, 0.2*u.um) + SpectralRegion(0.3*u.um, 0.4*u.um) +\
    ...      SpectralRegion(0.45*u.um, 0.6*u.um) + SpectralRegion(0.8*u.um, 0.9*u.um) +\
    ...      SpectralRegion(1.0*u.um, 1.2*u.um) + SpectralRegion(1.3*u.um, 1.5*u.um)

    >>> for s in sr:
    ...     print(s.lower) #doctest:+SKIP
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

    >>> sr = SpectralRegion(0.15*u.um, 0.2*u.um) + SpectralRegion(0.3*u.um, 0.4*u.um) +\
    ...      SpectralRegion(0.45*u.um, 0.6*u.um) + SpectralRegion(0.8*u.um, 0.9*u.um) +\
    ...      SpectralRegion(1.0*u.um, 1.2*u.um) + SpectralRegion(1.3*u.um, 1.5*u.um)

    >>> sr_inverted = sr.invert(0.05*u.um, 3*u.um)
    >>> sr_inverted #doctest:+SKIP
    Spectral Region, 7 sub-regions:
    (0.05 um, 0.15 um)   (0.2 um, 0.3 um)     (0.4 um, 0.45 um)
    (0.6 um, 0.8 um)     (0.9 um, 1.0 um)     (1.2 um, 1.3 um)
    (1.5 um, 3.0 um)

Region Extraction
-----------------

Given a `~specutils.SpectralRegion`, one can extract a sub-spectrum
from a `~specutils.Spectrum1D` object. If the `~specutils.SpectralRegion`
has multiple sub-regions then a list of `~specutils.Spectrum1D` objects will
be returned.

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

Reference/API
-------------

.. automodapi:: specutils
    :no-main-docstr:
    :no-heading:
    :no-inheritance-diagram:

    :skip: test
    :skip: Spectrum1D
    :skip: SpectrumCollection
    :skip: UnsupportedPythonError
    :skip: SpectralAxis
