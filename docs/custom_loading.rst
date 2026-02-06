*************************************************
Loading and Defining Custom Spectral File Formats
*************************************************

Loading From a File
-------------------

Specutils leverages the astropy io registry to provide an interface for conveniently
loading data from files. To create a custom loader, the user must define it in
a separate python file and place the file in their ``~/.specutils`` directory.

Loading from a FITS File
------------------------
A spectra with a *Linear Wavelength Solution* can be read using the ``read``
method of the :class:`~specutils.Spectrum` class to parse the file name and
format


.. code-block:: python

  import os
  from specutils import Spectrum

  file_path = os.path.join('path/to/folder', 'file_with_1d_wcs.fits')

  spec = Spectrum.read(file_path, format='wcs1d-fits')


This will create a :class:`~specutils.Spectrum` object that you can manipulate later.

For instance, you could plot the spectrum.

.. code-block:: python

  import matplotlib.pyplot as plt

  plt.title('FITS file with 1D WCS')
  plt.xlabel('Wavelength (Angstrom)')
  plt.ylabel('Flux (erg/cm2/s/A)')
  plt.plot(spec.wavelength, spec.flux)
  plt.show()


.. image:: img/read_1d.png


Creating a Custom Loader
------------------------

Defining a custom loader consists of importing the
`~specutils.io.registers.data_loader` decorator from specutils and attaching
it to a function that knows how to parse the user's data.  The return object
of this function must be an instance of one of the spectral classes
(:class:`~specutils.Spectrum`, :class:`~specutils.SpectrumCollection`,
:class:`~specutils.SpectrumList`).

Optionally, the user may define an identifier function. This function acts to
ensure that the data file being loaded is compatible with the loader function.

.. code-block:: python

    # ~/.specutils/my_custom_loader.py
    import os

    from astropy.io import fits
    from astropy.nddata import StdDevUncertainty
    from astropy.table import Table
    from astropy.units import Unit
    from astropy.wcs import WCS

    from specutils.io.registers import data_loader
    from specutils import Spectrum


    # Define an optional identifier. If made specific enough, this circumvents the
    # need to add ``format="my-format"`` in the ``Spectrum.read`` call.
    def identify_generic_fits(origin, *args, **kwargs):
        return (isinstance(args[0], str) and
                os.path.splitext(args[0].lower())[1] == '.fits')


    @data_loader("my-format", identifier=identify_generic_fits,
                 extensions=['fits'])
    def generic_fits(file_name, **kwargs):
        with fits.open(file_name, **kwargs) as hdulist:
            header = hdulist[0].header

            tab = Table.read(file_name)

            meta = {'header': header}
            wcs = WCS(hdulist[0].header)
            uncertainty = StdDevUncertainty(tab["err"])
            data = tab["flux"] * Unit("Jy")

        return Spectrum(flux=data, wcs=wcs, uncertainty=uncertainty, meta=meta)


An ``extensions`` keyword can be provided. This allows for basic filename
extension matching in the case that the ``identifier`` function is not
provided.

It is possible to query the registry to return the list of loaders associated
with a particular extension.

.. code-block:: python

    from specutils.io import get_loaders_by_extension

    loaders = get_loaders_by_extension('fits')

The returned list contains the format labels that can be fed into the ``format``
keyword argument of the ``Spectrum.read`` method.

After placing this python file in the user's ``~/.specutils`` directory, it
can be utilized by referencing its name in the ``read`` method of the
:class:`~specutils.Spectrum` class

.. code-block:: python

    from specutils import Spectrum

    spec = Spectrum.read("path/to/data", format="my-format")

.. _multiple_spectra:

Loading Multiple Spectra
^^^^^^^^^^^^^^^^^^^^^^^^

It is possible to create a loader that reads multiple spectra from the same
file. For the general case where none of the spectra are assumed to be the same
length, the loader should return a `~specutils.SpectrumList`. Consider the
custom JWST data loader as an example:

.. code-block:: python

    from specutils import Spectrum, SpectrumList
    from specutils.io.registers import data_loader

    @data_loader(
        "JWST x1d multi", identifier=identify_jwst_x1d_multi_fits,
        dtype=SpectrumList, extensions=['fits'], priority=10,
    )
    def jwst_x1d_multi_loader(file_obj, **kwargs):
        """Loader for JWST x1d 1-D spectral data in FITS format"""
        return _jwst_spec1d_loader(file_obj, extname='EXTRACT1D', **kwargs)

    def _jwst_spec1d_loader(file_obj, extname='EXTRACT1D', flux_col=None, **kwargs):
        """Implementation of loader for JWST x1d 1-D spectral data in FITS format"""

        if extname not in ['COMBINE1D', 'EXTRACT1D']:
            raise ValueError('Incorrect extname given for 1d spectral data.')

        spectra = []
        with read_fileobj_or_hdulist(file_obj, memmap=False, **kwargs) as hdulist:

            primary_header = hdulist["PRIMARY"].header

            for hdu in hdulist:
                # Read only the BinaryTableHDUs named COMBINE1D/EXTRACT1D and SCI
                if hdu.name != extname:
                    continue

                header = hdu.header

                # Correct some known bad unit strings before reading the table
                bad_units = {"(MJy/sr)^2": "MJy2 sr-2"}
                for c in hdu.columns:
                    if c.unit in bad_units:
                        c.unit = bad_units[c.unit]

                data = QTable.read(hdu)

                if data[0]['WAVELENGTH'].shape != ():
                    # In this case we have multiple spectra packed into a single extension, one target
                    # per row of the table
                    for row in data:
                        if hasattr(row['WAVELENGTH'], 'mask') and np.all(row['WAVELENGTH'].mask):
                            # If everything is masked out we don't bother to read it in at all
                            continue
                        srctype = row['SOURCE_TYPE']
                        spec = _jwst_spectrum_from_table(row, header, primary_header, flux_col, srctype)
                        spectra.append(spec)
                else:
                    # Otherwise the whole table is defining a single spectrum
                    spec = _jwst_spectrum_from_table(data, header, primary_header, flux_col)
                    spectra.append(spec)

        return SpectrumList(spectra)

Note that by default, any loader that uses ``dtype=Spectrum`` will also
automatically add a reader for `~specutils.SpectrumList`. This enables user
code to call `specutils.SpectrumList.read <astropy.nddata.NDIOMixin.read>` in
all cases if it can't make assumptions about whether a loader returns one or
many `~specutils.Spectrum` objects. This method is available since
`~specutils.SpectrumList` makes use of the Astropy IO registry (see
`astropy.io.registry.read`).

Lazy Loading
^^^^^^^^^^^^

By default, `~specutils.SpectrumList` data loaders will load all spectra eagerly.  Loaders optionally support
lazy loading so that individual spectra are only loaded into memory when accessed.  This can be useful for lists with a
large number of spectra.

Implementation
~~~~~~~~~~~~~~
Lazy loading is opt-in per data loader. To implement lazy loading for a given data loader, define a custom function to be
passed into the ``lazy_loader`` argument of the ``@data_loader`` decorator. The loader function should return a ``SpectrumList`` built using
the :meth:`specutils.SpectrumList.from_lazy` class method.  The function should:

* Determine the total number of spectra.
* Define an index-based loader function that returns a single ``Spectrum``.
* Return the resulting ``SpectrumList``.

See the following example for the Roman 1d spectra asdf data loader.

.. code-block:: python

    def _lazy_loader(file_obj, **kwargs):
        """Lazy loader for Roman spectra"""
        # read in the input file
        with read_fileobj_or_asdftree(file_obj, **kwargs) as af:
            roman = af["roman"]
            # get the roman spectral source ids
            sources = list(roman["data"].keys())

        def _loader(i: int) -> Spectrum:
            """Function to load a single spectra from the input file given a list index"""
            # select the proper source
            source = sources[i]
            with read_fileobj_or_asdftree(file_obj, **kwargs) as af2:
                roman2 = af2["roman"]
                # load a single Spectrum
                return _load_roman_spectrum(roman2, source)

        # create the lazy SpectrumList, pass in the number of spectra and the individual spectrum loader
        sl = SpectrumList.from_lazy(length=len(sources), loader=_loader)
        return sl


    @data_loader(
        "Roman 1d combined",
        identifier=identify_1d_combined,  # standard function for format identification
        dtype=SpectrumList,
        extensions=["asdf"],
        priority=10,
        force=True,
        lazy_loader=_lazy_loader,  # function to handle lazy loading
    )
    def roman_1d_combined_list(file_obj, **kwargs):
        """Load all Roman 1d combined extracted spectra"""
        # standard eager loading of all spectra
        spectra = SpectrumList()
        with read_fileobj_or_asdftree(file_obj, **kwargs) as af:
            roman = af["roman"]
            meta = roman["meta"]
            # load the spectra
            for source in roman["data"]:
                # load single spectrum
                spectrum = _load_roman_spectrum(roman, source)
                spectra.append(spectrum)

        return spectra


    def _load_roman_spectrum(roman: dict, source: str) -> Spectrum:
        """Load a single Roman spectrum"""
        meta = copy.deepcopy(roman["meta"])
        meta['source_id'] = source
        data = roman["data"][source] if source else roman["data"]
        flux = data['flux'] * u.Unit(meta["unit_flux"])
        flux_err  = StdDevUncertainty(data['flux_error'])
        wavelength = data['wl'] * u.Unit(meta["unit_wl"])
        return Spectrum(spectral_axis=wavelength, flux=flux, uncertainty=flux_err, meta=meta)

Usage
~~~~~

Once implmemented, lazy loading can be activated by passing ``lazy_load=True`` to ``SpectrumList.read``.
This creates a list of placeholder objects of length equal to the number of spectra loaded into the list.
The ``repr`` indicates a lazy list with how many spectra are currently loaded into memory

.. code-block:: python

    # example file with 6 spectral sources
    speclist = SpectrumList.read("/path/to/roman.asdf", format="Roman 1d combined", lazy_load=True)

    speclist
    lazy list: 0 items loaded; access an index to load a spectrum:
    [<object object at 0x127280500>, <object object at 0x127280500>,
    <object object at 0x127280500>, <object object at 0x127280500>,
    <object object at 0x127280500>, <object object at 0x127280500>]

    # inspect the lazy list
    len(speclist)
    6

    # verify it is lazy
    speclist.is_lazy
    True

    # check how many are loaded
    speclist.n_loaded
    0

Accessing a list item will lazily load the corresponding spectrum into memory.

.. code-block:: python

    # access the first spectrum
    speclist[0]
    <Spectrum(flux=[nan ... nan] W / (nm m2) (shape=(275,), mean=0.00000 W / (nm m2)); spectral_axis=<SpectralAxis [ 750.          752.47542943  754.95902919 ... 1837.848077   1843.91402794
    1850.        ] nm> (length=275); uncertainty=StdDevUncertainty)>

    # check the repr again
    speclist
    lazy list: 1 items loaded; access an index to load a spectrum:
    [<Spectrum(flux=[nan ... nan] W / (nm m2) (shape=(275,), mean=0.00000 W / (nm m2)); spectral_axis=<SpectralAxis [ 750.          752.47542943  754.95902919 ... 1837.848077   1843.91402794
    1850.        ] nm> (length=275); uncertainty=StdDevUncertainty)>,
    <object object at 0x127280500>, <object object at 0x127280500>, <object object at 0x127280500>,
    <object object at 0x127280500>, <object object at 0x127280500>]

    # check how many are loaded
    speclist.n_loaded
    1

**Optional Labels**
You can optionally pass a list of labels to use as placeholder values in the lazy list repr, instead of
the default pointer ``<object>``.  This can be done by passing a list of strings to the ``labels`` argument
of ``SpectrumList.from_lazy`` in the lazy loader function.

.. code-block:: python

    def _lazy_load_roman(file_obj, **kwargs):
        """Lazy loader for SpectrumList"""

        with read_fileobj_or_asdftree(file_obj, **kwargs) as af:
            roman = af["roman"]
            # create a list of roman source ids
            sources = list(roman["data"].keys())

        def _loader(i: int) -> Spectrum:
            ...

        # pass the source ids as placeholder labels
        sl = SpectrumList.from_lazy(
            length=len(sources), loader=_loader, labels=sources
        )
        return sl

Loading the lazy list with display these labels instead:

.. code-block:: python

    speclist = SpectrumList.read("/path/to/roman.asdf", format="Roman 1d combined", lazy_load=True)

    speclist
    lazy list: 0 items loaded; access an index to load a spectrum:
    ['402849', '403613', '403686', '404935', '404979', '414981']


.. note::

    Lazy loaders can be outfitted to any existing data loader.  See the example data loader for loading
    ``SDSS-V spec`` formatted FITS files.


Alternate List Indexing
^^^^^^^^^^^^^^^^^^^^^^^

Alternate ID labels allow string indexing of a ``SpectrumList``.  This is useful for long lists of
spectra that can be more easily identified by a name or ID rather than a list index.  Alternate IDs
are optional, and can be added to any ``SpectrumList`` data loader with the :meth:`specutils.SpectrumList.set_id_map`
class method.  This method accepts a dictionary mapping of string labels to list indices.

For example,

.. code-block:: python

    from specutils import SpectrumList

    # instantiate a SpectrumList
    ss = SpectrumList(['a', 'b', 'c'])
    ss
    ['a', 'b', 'c']

    # set an alternate id indexing
    ss.set_id_map({'spec1': 0, 'spec2': 1, 'spec3': 2})

    # access an item with a list index
    ss[0]
    'a'

    # access an item with an alternate id
    ss['spec1']
    'a'

The example Roman data loaders uses a string target source id as alternate ids.

.. code-block:: python

    def _load_roman_multisource(file_obj, **kwargs):
        """Load all Roman spectra into a SpectrumList"""

        spectra = SpectrumList()
        with read_fileobj_or_asdftree(file_obj, **kwargs) as af:
            roman = af["roman"]
            meta = roman["meta"]
            sources = list(roman['data'].keys())

            # set the alternate ids to roman source ids
            source_idx_map = dict(zip(sources, range(len(sources))))
            spectra.set_id_map(source_idx_map)

            # load the spectra
            for source in roman["data"]:
                spectrum = _load_roman_spectrum(roman, source)
                spectra.append(spectrum)

        return spectra

.. _custom_writer:

Creating a Custom Writer
------------------------

Similar to creating a custom loader, a custom data writer may also be defined.
This again will be done in a separate python file and placed in the user's
``~/.specutils`` directory to be loaded into the astropy io registry.

.. code-block:: python

    # ~/.spectacle/my_writer.py
    from astropy.table import Table
    from specutils.io.registers import custom_writer


    @custom_writer("fits-writer")
    def generic_fits(spectrum, file_name, **kwargs):
        flux = spectrum.flux.value
        disp = spectrum.spectral_axis.value
        meta = spectrum.meta

        tab = Table([disp, flux], names=("spectral_axis", "flux"), meta=meta)

        tab.write(file_name, format="fits")

The custom writer can be used by passing the name of the custom writer to the
``format`` argument of the ``write`` method on the
:class:`~specutils.Spectrum`.

.. code-block:: python

    spec = Spectrum(flux=np.random.sample(100) * u.Jy,
                      spectral_axis=np.arange(100) * u.AA)

    spec.write("my_output.fits", format="fits-writer")


Reference/API
-------------
.. automodapi:: specutils.io.registers
    :no-heading:
