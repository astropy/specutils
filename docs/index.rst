
.. the "raw" directive below is used to hide the title in favor of just the logo being visible
.. raw:: html

    <style media="screen" type="text/css">
      h1 { display:none; }
    </style>

***********************
Specutils Documentation
***********************

.. py:function:: specutils

.. image:: img/logo.png


The ``specutils`` package provides a basic interface for the loading,
manipulation, and common forms of analysis of spectroscopic data. These
generic data containers and accompanying modules will provide a toolbox that the
astronomical community can use to build more domain-specific packages. For more
details about the underlying principles, see
`APE13 <https://github.com/astropy/astropy-APEs/blob/master/APE13.rst>`_, the
guiding document for spectroscopic development in the Astropy Project.

.. note::
    While specutils is available for general use, the API is in an early enough
    development stage that some interfaces may change if user feedback and
    experience warrants it.

Getting started with `specutils`
================================

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
    >>> quantity_support()  # for getting units on the axes below # doctest: +ELLIPSIS
    <...>

    Now we load the dataset from it's canonical source:

    >>> f = fits.open('https://dr14.sdss.org/optical/spectrum/view/data/format=fits/spec=lite?plateid=1323&mjd=52797&fiberid=12')  # doctest: +IGNORE_OUTPUT
    >>> specdata = f[1].data # The spectrum is in the second HDU of this file.
    >>> f.close()

    Then we re-format this dataset into astropy quantities, and create a
    `~specutils.Spectrum1D` object:

    >>> from specutils import Spectrum1D
    >>> lamb = 10**specdata['loglam'] * u.AA
    >>> flux = specdata['flux'] * 10**-17 * u.Unit('erg cm-2 s-1 AA-1')
    >>> spec = Spectrum1D(spectral_axis=lamb, flux=flux)

    And we plot it:

    >>> lines = plt.step(spec.spectral_axis, spec.flux)

Now maybe you want the equivalent width of a spectral line.  That requires
normalizing by a continuum estimate:

.. testsetup::

    >>> fig = plt.figure()  # necessary because otherwise the doctests fail due to quantity_support and the flux units being different from the last figure

.. plot::
    :include-source:
    :align: center
    :context: close-figs

    >>> from specutils.fitting import fit_generic_continuum
    >>> cont_norm_spec = spec / fit_generic_continuum(spec)(spec.spectral_axis)
    >>> lines = plt.step(cont_norm_spec.wavelength, cont_norm_spec.flux)
    >>> plt.xlim(654*u.nm, 660*u.nm)  # doctest: +FLOAT_CMP
    (6540., 6600.)

But then you can apply a single function over the region of the spectrum
containing the line:

.. code-block:: python

    >>> from specutils import SpectralRegion
    >>> from specutils.analysis import equivalent_width
    >>> equivalent_width(cont_norm_spec, regions=SpectralRegion(6562*u.AA, 6575*u.AA))
    <Quantity -14.78092438 Angstrom>




While there are other tools and spectral representations detailed more below,
this gives a test of the sort of analysis `specutils` enables.


Using `specutils`
=================

For more details on usage of specutils, see the sections listed below.

.. toctree::
    :maxdepth: 2

    installation
    types_of_spectra
    spectrum1d
    spectrum_collection
    spectral_regions
    analysis
    fitting
    manipulation
    arithmetic
    custom_loading
    contributing


Get Involved
------------

Please see :doc:`contributing` for information on bug reporting and
contributing to the specutils project.

.. toctree::
   :maxdepth: 2

   contributing
