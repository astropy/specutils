0.6 (unreleased)
----------------

API Changes
^^^^^^^^^^^


New Features
^^^^^^^^^^^^

Bug Fixes
^^^^^^^^^


Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

0.5.3
-----

Bug Fixes
^^^^^^^^^

- Fix comparison of FITSWCS objects in arithmetic operations.

- Fix example documentation when run in python interpreter.


0.5.2 (2019-02-06)
----------------

Bug Fixes
^^^^^^^^^
- Bugfixes for astropy helpers, pep8 syntax checking, and plotting in docs [#416,#417,#419]

- All automatically generated ``SpectrumList`` loaders now have identifiers. [#440]

- ``SpectralRange.from_center`` parameters corrected after change to SpectralRange interface. [#433]

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^
- Improve explanation on creating spectrum continua. [#420]

- Wrap IO identifier functions to ensure they always return True or False and log any errors. [#404]


0.5.1 (2018-11-29)
------------------

Bug Fixes
^^^^^^^^^

- Fixed a bug in using spectral regions that have been inverted. [#403]

- Use the pytest-remotedata plugin to control tests that require access to remote data. [#401,#408]


0.5 (2018-11-21)
----------------

This was the first release of specutils executing the
[APE14](https://github.com/astropy/astropy-APEs/blob/master/APE14.rst)
plan (i.e. the "new" specutils) and therefore intended for broad use.
