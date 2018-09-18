
.. the "raw" directive below is used to hide the title in favor of just the logo being visible
.. raw:: html

    <style media="screen" type="text/css">
      h1 { display:none; }
    </style>

***********************
Specutils Documentation
***********************

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
    arithmetic
    smoothing
    custom_loading
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
