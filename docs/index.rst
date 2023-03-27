
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


Changes in version 2
====================

Specutils version 2 implemented a major change in that `~specutils.Spectrum1D`
no longer forces the spectral axis to be last for multi-dimensional data. This
was motivated by the desire for greater flexibility to allow for interoperability
with other packages that may wish to use ``specutils`` classes as the basis for
their own, and by the desire for consistency with the axis order that results
from a simple ``astropy.io.fits.read`` of a file. The legacy behavior can be
replicated by setting ``move_spectral_axis='last'`` when creating a new
`~specutils.Spectrum1D` object.

For a summary of other changes in version 2, please see the
`release notes <https://github.com/astropy/specutils/releases>`_.


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
    <Quantity -14.7396 Angstrom>


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
