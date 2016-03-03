.. _doc_custom_loaders:

Custom Loaders
==============

SpecViz utilizes
`Astropy I/O registry <http://docs.astropy.org/en/latest/io/registry.html>`_
and `YAML data serialization language <http://yaml.org/>`_  to enable flexible
support for a variety of data formats both in FITS and ASCII.

When ``specviz`` is called with a filename as argument
(see :ref:`doc_launching`), the default formats below are assumed based on file
extension:

* ``.fits`` or ``.mits`` -- :ref:`doc_ref_fits_loader`
* ``.txt`` or ``.dat`` -- :ref:`doc_def_ascii_loader`

By examining the
`YAML definitions <https://github.com/spacetelescope/specviz/tree/master/specviz/interfaces/default_loaders>`_
in the following sub-sections and their associated
`example data files <https://github.com/spacetelescope/specviz/tree/master/specviz/data>`_,
you will be able to create your own YAML file to define most custom data formats
(see :ref:`doc_create_custom_loader`).
In addition, there are also other custom loaders that come with SpecViz that
follow the same rules but are modified to load JWST DADF test data, which you
can also use as a reference.

While using YAML is very flexible, it is also very sensitive to slight changes
in your file format. For instance, if you have two files from the same
instrument but processed differently (say, one was extracted using IRAF and
another one using your own IDL program), they might have different formats
and will need separate YAML files.


.. _doc_ref_fits_loader:

Generic FITS Loader
-------------------

::

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

This is
`generic.yaml <https://github.com/spacetelescope/specviz/blob/master/specviz/interfaces/default_loaders/generic.yaml>`_,
which is a built-in YAML definition for a
`"generic" FITS file <https://github.com/spacetelescope/specviz/blob/master/specviz/data/generic_spectra.fits>`_.

::

  --- !CustomLoader

The first line states that this YAML file defines our custom loader. This is
always the same no matter what kind of format you are defining.

::

  name: Generic Fits
  extension: [fits, mits]

These two lines define the format name and accepted extensions, respectively.
In SpecViz GUI, this will translate to "Generic Fits (\*.fits \*.mits)" in the
file type drop-down menu.

::

  wcs:
    hdu: 0

This instructs the loader to look for WCS information in the ``PRIMARY``
(Extension 0) header. SpecViz also uses this header for other look-ups, e.g.,
flux unit from ``BUNIT`` (see below). Therefore, even if there are no WCS
information in your file, you always define this block and point the HDU value
to the ``PRIMARY`` header. If WCS information are available and supported by
:ref:`Astropy WCS <astropy:astropy-wcs>`, they will be used to establish
dispersion values and unit.

In the absence of WCS or the presence of explicit dispersion column, an
additional ``dispersion:`` block (not shown) can be defined similarly as
``data:`` (see below). Unlike flux, ``TNULL`` masking is ignored. If its
column does not have ``TUNIT`` and ``unit: 'unitname'`` is defined, SpecViz
will fall back to WCS unit. If all these unit look-ups failed, it defaults to
null unit. In the presence of ``unit:`` definition, it overrides both ``TUNIT``
and WCS.

::

  data:
    hdu: 1
    col: 0

This instructs the loader to look for flux values (data) in Extension 1, the
first column (column index starts from 0). If the column complies to FITS
standards (see `Astropy FITS Table <http://docs.astropy.org/en/stable/io/fits/usage/table.html>`_),
flux unit (inferred from ``TUNIT``) and data mask (inferred from ``TNULL``) are
also extracted from the same column.

If ``TUNIT`` is not defined, loader will look for ``unit: 'unitname'``
definition within this block, the unit name must be one that is accepted by
:ref:`Astropy Units <astropy:astropy-units>` (case sensitive). If that is
undefined as well, flux unit is extracted from ``BUNIT`` keyword in the same
header that contains WCS information (see ``wcs:`` block for details).
If all unit look-ups failed, flux unit is assumed to be
:math:`\textnormal{erg} \; \AA^{-1} \; \textnormal{cm}^{-2} \; \textnormal{s}^{-1}`.

Note that ``TNULL`` is not the same as DQ arrays, which can be similarly defined
with ``mask:`` block (not shown). If both ``TNULL`` and DQ are defined, the
masks will be combined.

::

  uncertainty:
    hdu: 1
    col: 1
    type: 'std'

This instructs the loader to look for flux uncertainty values in Extension 1,
the second column. Uncertainty type ``'std'`` states that the values are
standard deviation (as opposed to inverse variance, ``'ivar'``). Unlike flux
data, its ``TNULL`` masking is ignored and ``unit:`` tag is not supported.
If ``TUNIT`` is present, loader will attempt to convert the values to flux unit
first. Otherwise, its unit is assumed to be the same as flux unit.
If inverse variance is given, square-root is applied to the inversed values
before being converted to `~astropy.nddata.StdDevUncertainty`.

::

  meta:
    author: Nicholas Earl

The ``meta:`` block can contain any metadata tags you wish to include. They do
not affect how SpecViz works. In this example, the ``author:`` tag identifies
Nicholas Earl as the origin author of this YAML file.


.. _doc_def_ascii_loader:

ASCII Loader
------------

::

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

This is
`ascii.yaml <https://github.com/spacetelescope/specviz/blob/master/specviz/interfaces/default_loaders/ascii.yaml>`_,
which is a built-in YAML definition for a
`"generic" ASCII file <https://github.com/spacetelescope/specviz/blob/master/specviz/data/generic_spectra.txt>`_.

::

  --- !CustomLoader

The first line states that this YAML file defines our custom loader. This is
always the same no matter what kind of format you are defining.

::

  name: ASCII
  extension: [txt, dat]

These two lines define the format name and accepted extensions, respectively.
In SpecViz GUI, this will translate to "ASCII (\*.txt \*.dat)" in the
file type drop-down menu. All ASCII files must comply to
:ref:`Astropy ASCII Table <astropy:io-ascii>` standards.

Any header comments with ``KEY = VALUE`` format will be extracted as header
metadata information (currently not used by SpecViz).

::

  dispersion:
    col: 0
    unit: 'Angstrom'

Unlike :ref:`doc_ref_fits_loader`, ASCII table does not contain WCS. Therefore,
the ``dispersion:`` block is necessary to define the actual dispersion
(e.g., wavelength) values. This instructs the loader to look for dispersion
values in the first column (column index starts from 0).
Its unit, if not defined in the table itself
(e.g., via `IPAC table format <http://irsa.ipac.caltech.edu/applications/DDGEN/Doc/ipac_tbl.html>`_),
will be taken from the ``unit:`` tag. The given unit name must be one that is
accepted by :ref:`Astropy Units <astropy:astropy-units>` (case sensitive).
If unit is defined in both table and tag, the latter is ignore.
If unit is not defined anywhere, it defaults to Angstrom.

::

  data:
    col: 1
    unit: 'erg / (Angstrom cm2 s)'

This instructs the loader to look for flux values (data) in the second column.
Flux unit handling is similar to dispersion unit (see above), except that the
default unit would be
:math:`\textnormal{erg} \; \AA^{-1} \; \textnormal{cm}^{-2} \; \textnormal{s}^{-1}`
if undefined.

If there is an associated DQ column, it can be extracted in a similar fashion
using a ``mask:`` block specifying the column index (unit is not applicable).
Like :ref:`doc_ref_fits_loader`, zero mask values signify good data.

::

  uncertainty:
    col: 2
    type: 'std'

This instructs the loader to look for flux uncertainty values in the third
column. Uncertainty type ``'std'`` states that the values are standard deviation
(as opposed to inverse variance, ``'ivar'``). Its unit must be the same as flux
unit. If inverse variance is given, square-root is applied to the inversed
values before being converted to `~astropy.nddata.StdDevUncertainty`.

::

  meta:
    author: STScI

The ``meta:`` block can contain any metadata tags you wish to include. They do
not affect how SpecViz works. In this example, the ``author:`` tag identifies
STScI as the origin author of this YAML file.


.. _doc_create_custom_loader:

Creating a Custom Loader
------------------------

In addition to loaders that come pre-packaged with the software, SpecViz also
looks for custom loaders that you created and saved in your ``~/.specviz``
directory, which can be created with the following Unix command::

    $ mkdir ~/.specviz

To create your own loader, you can use either :ref:`doc_ref_fits_loader` or
:ref:`doc_def_ascii_loader` as a template. Your YAML file can have any name of
your choosing but must end with a ``.yaml`` extension.

In this section, we use a `MOSFIRE <http://www2.keck.hawaii.edu/inst/mosfire/>`_
spectrum named
`spec1d.gds1312_H0.003.emp26177.fits <https://github.com/spacetelescope/specviz/tree/master/specviz/data/spec1d.gds1312_H0.003.emp26177.fits>`_
as an example of a custom FITS format for which we must create our own custom
YAML definition file from the :ref:`doc_ref_fits_loader` template.

First, we inspect the file format that we have, as follow.

.. code-block:: python

    >>> from astropy.io import fits
    >>> pf = fits.open('spec1d.gds1312_H0.003.emp26177.fits')
    >>> pf.info()
    Filename: spec1d.gds1312_H0.003.emp26177.fits
    No.    Name         Type      Cards   Dimensions   Format
    0    PRIMARY     PrimaryHDU       4   ()
    1                BinTableHDU     23   1R x 3C      [2287E, 2287E, 2287E]

This opens the FITS file and prints out the overall file structure. From this,
it is obvious that the table is in Extension 1.

    >>> pf[0].header
    SIMPLE  =                    T /Dummy Created by MWRFITS v1.4a
    BITPIX  =                    8 /Dummy primary header created by MWRFITS
    NAXIS   =                    0 /No data is associated with this header
    EXTEND  =                    T /Extensions may (will!) be present
    >>> pf[1].header
    XTENSION= 'BINTABLE'           /Binary table written by MWRFITS v1.4a
    BITPIX  =                    8 /Required value
    ...
    TFORM3  = '2287E   '           /

This prints all the headers and we find no WCS information in either of the
extensions.

    >>> from astropy.table import Table
    >>> tab = Table.read(pf[1], format='fits')
    >>> print(tab)
    FLUX [2287]      LAMBDA [2287]      IVAR [2287]
    ---------------- ------------------ ---------------
    0.0 .. 0.0989667 14500.0 .. 18223.8 1e-06 .. 422.54

This shows that there are three columns in the table in Extension 1, namely
flux, wavelength, and inverse variance. The table has 2287 rows. Knowing the
wavelength regime that the instrument is sensitive to and looking at the
wavelength values, we can safely assume that the wavelength unit is Angstrom.

    >>> from astropy import units as u
    >>> u.electron / u.s / u.pix
    Unit("electron / (pix s)")

However, the flux unit is not defined anywhere and cannot be easily inferred.
So, let's just say that we already know the unit to be electrons/s/pix. The code
above shows us how Astropy can ingest the flux unit that we want.

::

  --- !CustomLoader
  name: Keck/MOSFIRE Fits
  extension: [fits, mits]
  wcs:
    hdu: 0
  dispersion:
    hdu: 1
    col: 1
    unit: 'Angstrom'
  data:
    hdu: 1
    col: 0
    unit: 'electron / (pix s)'
  uncertainty:
    hdu: 1
    col: 2
    type: 'ivar'
  meta:
    author: STScI

Now that we have the format figured out, it is time to write our own YAML file
for it. We will name it
`keck_mosfire.yaml <https://github.com/spacetelescope/specviz/blob/master/specviz/interfaces/default_loaders/keck_mosfire.yaml>`_.

::

  --- !CustomLoader

The first line states that this YAML file defines our custom loader.

::

  name: Keck/MOSFIRE Fits
  extension: [fits, mits]

These two lines define the format name and accepted extensions, respectively.
We will keep the extensions from our FITS template but change the name to
identify our new format. In SpecViz GUI, this will translate to
"Keck/MOSFIRE Fits (\*.fits \*.mits)" in the file type drop-down menu.

::

  wcs:
    hdu: 0

We do not have WCS nor ``BUNIT`` defined, so we will simply leave this the same
as our template.

::

  dispersion:
    hdu: 1
    col: 1
    unit: 'Angstrom'

This instructs the loader to look for our wavelength values in Extension 1, the
second column. We explicitly set its unit to Angstrom.

::

  data:
    hdu: 1
    col: 0
    unit: 'electron / (pix s)'

This instructs the loader to look for our flux values in Extension 1, the first
column, like the template. However, we also explicitly set its unit to
electrons/s/pix by providing the appropriate Astropy unit name.

::

  uncertainty:
    hdu: 1
    col: 2
    type: 'ivar'

This instructs the loader to look for flux uncertainty values in Extension 1,
the third column. Unlike the template, we define it as inverse variance.

::

  meta:
    author: STScI

Since this does not affect how SpecViz works, we do the lazy thing here by
leaving it the same as our template.

Once you are done writing your YAML file, be sure to save it in ``~/.specviz``.
Next, start SpecViz as usual. Now, in its open file dialog, you will see
your new format listed in the file-type drop-down menu.
