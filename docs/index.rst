
.. the "raw" directive below is used to hide the title in favor of just the logo being visible
.. raw:: html

    <style media="screen" type="text/css">
      h1 { display:none; }
    </style>

***********************
Specutils Documentation
***********************

.. _specutils:

.. image:: _static/logo.png


``specutils`` is a Python package for representing, loading,
manipulating, and analyzing astronomical spectroscopic data. The
generic data containers and accompanying modules provide a toolbox that the
astronomical community can use to build more domain-specific packages. For more
details about the underlying principles, see
`APE13 <https://github.com/astropy/astropy-APEs/blob/main/APE13.rst>`_, the
guiding document for spectroscopic development in the Astropy Project.

.. note::
    While specutils is available for general use, the API is in an early enough
    development stage that some interfaces may change if user feedback and
    experience warrants it.

Changes coming in version 2.0
=============================

Specutils 2.0 has been in development for some time and is nearly ready for release.
The major changes that will affect users are detailed here in an attempt to prepare
users for the transition.

The most visible change is that the `~specutils.Spectrum1D` class will be renamed
to ``Spectrum`` to reduce confusion about multi-dimensional flux arrays being supported.
The current class name will be deprecated in version 2.1; importing the old name will
work but raise a deprecation warning until then. Version 1.20 implemented a ``Spectrum``
class as a simple wrapper around `~specutils.Spectrum1D` so that you may update your
code to the new class name now and avoid deprecation warnings when 2.0 releases. Note
that the new keyword arguments ``move_spectral_axis`` and ``spectral_axis_index`` being
introduced in 2.0 will be ignored in 1.x if used when initializing the ``Spectrum`` class.

Single-dimensional flux use cases should be mostly unchanged in 2.0, with the exception
being that spectrum arithmetic will check that the spectral axis of both operands are
equal, rather than simply checking that they are the same length. Thus, you will need
to resample onto a common spectral axis if doing arithmetic on spectra with differing
spectral axes.

Specutils version 2 implements a major change in that ``Spectrum``
no longer forces the spectral axis to be last for multi-dimensional data. This
was motivated by the desire for greater flexibility to allow for interoperability
with other packages that may wish to use ``specutils`` classes as the basis for
their own, and by the desire for consistency with the axis order that results
from a simple ``astropy.io.fits.read`` of a file. The legacy behavior can be
replicated by setting ``move_spectral_axis='last'`` when creating a new
``Spectrum`` object. ``Spectrum`` will attempt to automatically
determine which flux axis corresponds to the spectral axis during initialization
based on the WCS (if provided) or the shape of the flux and spectral axis arrays,
but if the spectral axis index is unable to be automatically determined you will
need to specify which flux array axis is the dispersion axis with the
``spectral_axis_index`` keyword. Note that since the ``spectral_axis`` can specify
either bin edges or bin centers, a flux array of shape ``(10, 11)`` with spectral axis
of length 10 or 11 would be ambigious. In this case you could initialize a
``Spectrum`` with ``bin_specification`` set to either "edges" or "centers"
to break the degeneracy.

An additional change for multi-dimensional spectra is that previously, initializing
such a ``Spectrum`` with a  ``spectral_axis`` specified, but no WCS, would
create a ``Spectrum`` instance with a one-dimensional GWCS that was essentially
a lookup table with the spectral axis values. In 2.0 this case will result in a GWCS with
dimensionality matching that of the flux array to facilitate use with downstream packages
that expect WCS dimensionality to match that of the data. The resulting spatial axes
transforms are simple pixel to pixel identity operations, since no actual spatial
coordinate information is available.

In addition to the changes to the generated GWCS, handling of input GWCS will also be
improved. This mostly manifests in the full GWCS (including spatial information) being
retained in the resulting ``Spectrum`` objects when reading, e.g., JWST spectral
cubes.

Getting started with :ref:`specutils <specutils>`
=================================================

As a basic example, consider an emission line galaxy spectrum from the
`SDSS <https://www.sdss.org/>`_.  We will use this as a proxy for a spectrum you
may have downloaded from some archive, or reduced from your own observations.

.. plot::
    :include-source:
    :align: center
    :context: close-figs

    We begin with some basic  imports:

    >>> from astropy.io import fits
    >>> from astropy import units as u
    >>> import numpy as np
    >>> from matplotlib import pyplot as plt
    >>> from astropy.visualization import quantity_support
    >>> quantity_support()  # for getting units on the axes below  # doctest: +IGNORE_OUTPUT

    Now we load the dataset from its canonical source:

    >>> filename = 'https://data.sdss.org/sas/dr16/sdss/spectro/redux/26/spectra/1323/spec-1323-52797-0012.fits'
    >>> # The spectrum is in the second HDU of this file.
    >>> with fits.open(filename) as f:  # doctest: +IGNORE_OUTPUT +REMOTE_DATA
    ...     specdata = f[1].data  # doctest: +REMOTE_DATA

    Then we re-format this dataset into astropy quantities, and create a
    `~specutils.Spectrum1D` object:

    >>> from specutils import Spectrum1D
    >>> lamb = 10**specdata['loglam'] * u.AA # doctest: +REMOTE_DATA
    >>> flux = specdata['flux'] * 10**-17 * u.Unit('erg cm-2 s-1 AA-1') # doctest: +REMOTE_DATA
    >>> spec = Spectrum1D(spectral_axis=lamb, flux=flux) # doctest: +REMOTE_DATA

    And we plot it:

    >>> f, ax = plt.subplots()  # doctest: +IGNORE_OUTPUT +REMOTE_DATA
    >>> ax.step(spec.spectral_axis, spec.flux) # doctest: +IGNORE_OUTPUT +REMOTE_DATA

.. testsetup::

    >>> fig = plt.figure()  # necessary because otherwise the doctests fail due to quantity_support and the flux units being different from the last figure

.. plot::
    :include-source:
    :align: center
    :context: close-figs

    Now maybe you want the equivalent width of a spectral line. That requires
    normalizing by a continuum estimate:

    >>> import warnings
    >>> from specutils.fitting import fit_generic_continuum
    >>> with warnings.catch_warnings():  # Ignore warnings
    ...     warnings.simplefilter('ignore')
    ...     cont_norm_spec = spec / fit_generic_continuum(spec)(spec.spectral_axis) # doctest: +REMOTE_DATA

    >>> f, ax = plt.subplots()  # doctest: +IGNORE_OUTPUT +REMOTE_DATA
    >>> ax.step(cont_norm_spec.wavelength, cont_norm_spec.flux)  # doctest: +IGNORE_OUTPUT +REMOTE_DATA
    >>> ax.set_xlim(654 * u.nm, 660 * u.nm)  # doctest: +IGNORE_OUTPUT +REMOTE_DATA

    But then you can apply a single function over the region of the spectrum
    containing the line:

    >>> from specutils import SpectralRegion
    >>> from specutils.analysis import equivalent_width
    >>> equivalent_width(cont_norm_spec, regions=SpectralRegion(6562 * u.AA, 6575 * u.AA)) # doctest: +REMOTE_DATA +FLOAT_CMP
    <Quantity -14.82013888 Angstrom>


While there are other tools and spectral representations detailed more below,
this gives a test of the sort of analysis :ref:`specutils <specutils>` enables.


Using :ref:`specutils <specutils>`
==================================

For more details on usage of specutils, see the sections listed below.

.. toctree::
    :maxdepth: 2

    installation
    types_of_spectra
    spectrum1d
    spectrum_collection
    spectral_cube
    spectral_regions
    analysis
    fitting
    manipulation
    arithmetic
    wcs_utils
    custom_loading
    identify

Get Involved - Developer Docs
-----------------------------

Please see :doc:`contributing` for information on bug reporting and
contributing to the specutils project.

.. toctree::
   :maxdepth: 2

   contributing
   releasing
