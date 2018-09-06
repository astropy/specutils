.. image:: img/logo.png

The specutils package provides a basic interface for the loading, manipulation,
and low-level analysis of spectroscopic data. The intention is for the
generic data containers and accompanying modules to provide a basis for the
astronomical community upon which more robust tooling can be built.

The basic data container of specutils is the :class:`~specutils.Spectrum1D`
class. It contains logic to handle multi-dimensional flux data, spectral axis
in various forms (wavelenth, frequency, energy, velocity, etc.), convenient and
unobtrusive wcs support, and uncertainty handling. This core container also
maintaines a few other features including unit support, custom data loading,
arithmetic operation support, and multi-dimensional flux data support.

.. note:: support for collections of spectra with distinct spectral axes
          information is slotted for a future release.


Using specutils
---------------

.. toctree::
   :maxdepth: 2

   contributing
   installation
   getting_started
   custom_loading
   basic_analysis
   arithmetic
   smoothing
   fitting

.. toctree::
    :maxdepth: 1

    high-level_API.rst


Get Involved
------------


Please see :doc:`contributing` for information on bug reporting and
contributing to the specutils project.
