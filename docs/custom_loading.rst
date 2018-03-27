**************
Custom Loading
**************

Loading From a File
-------------------

Specutils leverages the astropy io registry to provide an interface for conveniently
loading data from files. To create a custom loader, the user must define it in
a separate python file and place the file in their `~/.specutils` directory.


Creating a Custom Loader
------------------------

Defining a custom loader consists of importing the `data_loader` decorator from
specutils and attaching it to a function that knows how to parse the user's data.
The return object of this function must be a :class:`~specutils.spectra.Spectrum1D`
object.

Optionally, the user may define an identifier function. This function acts to
ensure that the data file being loaded is compatible with the loader function.

.. code-block:: python

    # ~/.specutils/my_custom_loader.py
    import os
    import six

    from astropy.io import fits
    from astropy.nddata import StdDevUncertainty
    from astropy.table import Table
    from astropy.units import Unit
    from astropy.wcs import WCS

    from specutils.io.registers import data_loader
    from specutils import Spectrum1D

    # Define an optional identifier. If made specific enough, this circumvents the
    # need to add `format="my-format"` in the `Spectrum1D.read` call.
    def identify_generic_fits(origin, *args, **kwargs):
        return (isinstance(args[0], six.string_types) and
                os.path.splitext(args[0].lower())[1] == '.fits')


    @data_loader("my-format", identifier=identify_generic_fits)
    def generic_fits(file_name, **kwargs):
        name = os.path.basename(file_name.rstrip(os.sep)).rsplit('.', 1)[0]

        with fits.open(file_name, **kwargs) as hdulist:
            header = hdulist[0].header

            tab = Table.read(file_name)

            meta = {'header': header}
            wcs = WCS(hdulist[0].header)
            uncertainty = StdDevUncertainty(tab["err"])
            data = tab["flux"] * Unit("Jy")

        return Spectrum1D(flux=data, wcs=wcs, uncertainty=uncertainty, meta=meta)


After placing this python file in the user's `~/.specutils` directory, it can
be utilized by referencing its name in the `read` method of the :class:`~specutils.spectra.Spectrum1D`
class

.. code-block:: python

    from specutils import Spectrum1D

    spec = Spectrum1D.read("path/to/data", format="my-format")


Creating a Custom Writer
------------------------

Similar to creating a custom loader, a custom data writer may also be defined.
This again will be done in a separate python file and placed in the user's
`~/.specutils` directory to be loaded into the astropy io registry.

.. code-block:: python

    # ~/.spectacle/my_writer.py
    from astropy.table import Table
    from specutils.io.registers import custom_writer


    @custom_writer("fits-writer")
    def generic_fits(spectrum, file_name, **kwargs):
        flux = spectrum.flux.value
        disp = spectrum.dispersion.value
        meta = spectrum.meta

        tab = Table([disp, flux], names=("dispersion", "flux"), meta=meta)

        tab.write(file_name, format="fits")

The custom writer can be used by passing the name of the custom writer to the
`format` argument of the `write` method on the
:class:`~specutils.spectra.Spectrum1D`.

.. code-block:: python

    spec = Spectrum1D(flux=np.random.sample(100) * u.Jy,
                      spectral_axis=np.arange(100) * u.AA)

    spec.write("my_output.fits", format="fits-writer")