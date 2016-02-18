.. _doc_custom_loaders:

Custom Loaders
==============

Pyfocal utilizes
`Astropy I/O registry <http://docs.astropy.org/en/latest/io/registry.html>`_
and `YAML data serialization language <http://yaml.org/>`_  to enable flexible
support for a variety of data formats both in FITS and ASCII.

For instance, this built-in YAML definition for Generic FITS below states that::

  --- !CustomLoader
  name: Generic Fits
  extension: [fits, mits]
  wcs:
    hdu: 0
  data:
    hdu: 1
    col: 0
  uncertainty:
    hdu: 1
    col: 1
    type: 'std'
  meta:
    author: Nicholas Earl

* Relevant file extensions are either ``.fits`` or ``.mits``.
* WCS information are in the ``PRIMARY`` (Extension 0) header. They will be
  used to establish dispersion values and unit.
* Flux values (data) are in Extension 1, the first column (column index starts
  from 0). (If flux unit not found in the table column definition, Pyfocal
  will attempt to get the unit from the ``BUNIT`` keyword in the same header
  that contains WCS information. If still not found, it is assumed to be
  :math:`\textnormal{erg} \; \AA^{-1} \; \textnormal{cm}^{-2} \; \textnormal{s}^{-1}`.)
* Flux uncertainties are also in Extension 1, the second column. The values are
  standard deviation (as opposed to variance, ``'ivar'``). If unit is present,
  it will be converted to flux unit. Otherwise, it is assumed to be the same
  as flux unit.
* The definition file was written by Nicholas Earl.

Meanwhile, this build-in YAML definition for ASCII below states that::

  --- !CustomLoader
  name: ASCII
  extension: [txt, dat]
  dispersion:
    col: 0
    unit: 'Angstrom'
  data:
    col: 1
    unit: 'erg / (Angstrom cm2 s)'
  uncertainty:
    col: 2
    type: 'std'
  meta:
    author: STScI

* Relevant file extensions are either ``.txt`` or ``.dat``.
* Dispersion (e.g., wavelength) values are in the first column (column index
  starts from 0). Its unit, if not explicitly defined (e.g., via IPAC header),
  will be set to Angstrom, which is also the default unit Pyfocal will use if
  no unit information is given.
* Flux values (data) are in the second column. Its unit, if not explicitly
  defined (e.g., via IPAC header), will be set to
  :math:`\textnormal{erg} \; \AA^{-1} \; \textnormal{cm}^{-2} \; \textnormal{s}^{-1}`,
  which is also the default unit Pyfocal will use if no unit information is
  given.
* Flux uncertainties are in the third column. The values are standard deviation
  (as opposed to variance). Its unit must be the same as flux values.
* The definition file was written by STScI.

To add support for a new file format (e.g., reading spectra from Extension 4
instead of Extension 1), user will only need to provide a new YAML definition
file without needing to modify any codes in Pyfocal.

.. note:: Support for user-defined YAML is planned for a future release.
