Writing Spectra to a FITS file
=================================

.. currentmodule:: specutils.io.write_fits

.. warning::
    FITS writer doesn't work well after spectral slicing yet. In addition, reading fits files and writing the
     same files again might not give the exact same files, although the intended spectra will remain the same.

This page explains how to write spectra to FITS format. For more information on reading spectra and FITs in general, go to
the :doc:`FITS reader documentation page<read_fits>`.

Currently the FITS writer only works with 1-D spectra objects. These objects include the linear spectra objects, or the
multispec spectra objects (which are stored as a list of spectra). `write_fits.write` method deciphers the type of object passed
 automatically and writes spectra to the given file in FITS format. An example to write a simple linear FITS file::

   >>> from specutils.io import read_fits, write_fits
   >>> myspec = read_fits.read_fits_spectrum1d('myfile.fits')
   >>> write_fits.write(myspec, 'nynewfile.fits')

Writing to a file in multispec FITS format is similar::

   >>> from specutils.io import read_fits, write_fits
   >>> myspec = read_fits.read_fits_spectrum1d('myfile.fits')
   >>> write_fits.write(myspec, 'mymultispec.fits')



Reading FITS "multispec" format WCS
------------------------------------

The multispec format holds multiple one-dimensional spectra in a single file. The current 1-D readers can be used to read
such a file. A list of spectra is returned whenever the input file is a multispec file.
Here is an example of reading a simple FITS multispec format::

    >>> from specutils.io import read_fits
    >>> read_fits.read_fits_spectrum1d('mymultispec.fits')

The multispec format supports various functions to map the pixel indices to dispersion. Currently the following multispec formats
are supported:

 * Linear Dispersion Function
 * Log-linear Dispersion Function
 * Legendre Polynomial Dispersion Function
 * Chebyshev Polynomial Dispersion Function
 * Linear Spline Dispersion Function
 * Cubic Spline Polynomial Dispersion Function


Reference/API
-------------

.. automodapi:: specutils.io.write_fits
    :no-inheritance-diagram:
    :skip: OrderedDict, Spectrum1D
