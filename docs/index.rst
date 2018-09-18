.. image:: img/logo.png

The ``specutils`` package provides a basic interface for the loading,
manipulation, and low-level analysis of spectroscopic data. The intention is for
the generic data containers and accompanying modules to provide a basis for the
astronomical community upon which additional more domain-specific tools can be
built. To see more details about the underlying principals, see
`APE13 <https://github.com/astropy/astropy-APEs/blob/master/APE13.rst>`_, the
guiding document for spectroscopic development in the Astropy Project.



Types of Spectral Representations
---------------------------------

The core principle of ``specutils`` is to work as a toolbox.  That is, it aims
to provide the pieces needed to build particular spectroscopic workflows
without imposing a specific *required* set of algorithms or approaches to
spectroscopic analysis.  To that end, it aims to represent several different
types of ways one might wish to represent sets of spectroscopic data in Python.
These objects contains logic to handle multi-dimensional flux data, spectral
axes in various forms (wavelenth, frequency, energy, velocity, etc.), convenient
and unobtrusive wcs support, and uncertainty handling. The core containers also
handle units, a framework for reading and writing from various file formats,
arithmetic operation support, and a variety of analysis and manipulation tools
that work on them.


The core data objects  are primarily distinguished by the different ways of
storing the flux and the spectral axis . These cases are detailed below, along
with their corresponding ``specutils`` representations:

1. A 1D flux of length ``n``, and a matched spectral axis (which may be
   tabulated as an array, or encoded in a WCS). This is what typically is in
   mind when one speaks of "a single spectrum", and therefore the analysis tools
   are general couched as applying to this case. In ``specutils`` this is
   represented by the `~specutils.Spectrum1D` object with a 1-dimensional
   ``flux``.
2. A set of fluxes that can be represented in an array-like form of shape
   ``n x m (x ...)``,  with a spectral axis strictly of length ``n`` (and a
   matched WCS). In ``specutils`` this is represented by the
   `~specutils.Spectrum1D` object where ``len(flux.shape) > 1`` . In this sense
   the "1D" refers to the spectral axis, *not* the flux.
3. A set of fluxes  of shape ``n x m (x ...)``, and a set of spectral axes that
   are the same shape. This is distinguished from the above cases because there
   are as many spectral axes as there are spectra.  In this sense it is a
   collection of spectra, so can be thought of as a collection of
   `~specutils.Spectrum1D` objects.  But because it is often more performant to
   store the collection together as just one set of flux and spectral axis
   arrays, this case is represented by a separate object in ``specutils``:
   `~specutils.SpectrumCollection`.
4. An arbitrary collection of fluxes that are not all the same spectral length
   even in the spectral axis direction.  That is, case 3, but "ragged" in the
   sense that not all the spectra are length ``n``.  Because there is no
   performance benefit to be gained from using arrays (because the flux array is
   not rectangular), this case does not have a specific representation in
   ``specutils``.  Instead, this case should be dealt with by making lists (or
   numpy object-arrays) of `~specutils.Spectrum1D` objects, and iterating over
   them.

In all of these cases, the objects have additional attributes (e.g.
uncertainties), along with other metadata.  But the list above is exhaustive
under the assumption that the additional attributes have matched shape to either
flux or spectral axis (or some combination of the two).



Using specutils
---------------

.. toctree::
    :maxdepth: 2

    installation
    getting_started
    spectrum_collection
    spectral_regions
    custom_loading
    basic_analysis
    arithmetic
    smoothing
    fitting
    contributing

.. toctree::
    :maxdepth: 1

    high-level_API.rst


Get Involved
------------

Please see :doc:`contributing` for information on bug reporting and
contributing to the specutils project.

.. toctree::
   :maxdepth: 2

   contributing
