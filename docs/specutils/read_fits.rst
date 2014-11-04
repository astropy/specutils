Reading Spectra from a FITS file
=================================

.. currentmodule:: specutils.io.read_fits


FITS is one of the most common ways to store spectra. In most cases the flux values are stored as the pixel values
in a 1-, 2- or 3-D arrays. The mapping of axes indices to physically meaningful values (e.g. to wavelength) is encoded
in so FITS keywords (see :doc:`WCS <specwcs>`). This document describes how to extract this information and read it in to the various spectral
classes. A lot of information about the keyword header storage is taken from
`The IRAF/NOAO Spectral World Coordinate Systems <http://iraf.net/irafdocs/specwcs.php>`_. As the FITS keyword headers
have a complex way of storing information (often redundantly), we have created a FITS WCS parser (`~specutils.io.read_fits.FITSWCSSpectrum`)
that initializes with a FITS header and can extract/validate the information stored in the various keyword headers
making them available in a simple API. All of the FITS readers use this parser object.

.. note::
    FITS keywords can encode two mappings from physical to logical to dispersion. A general approach to logical coordinate
    systems is in development. Currently, the FITS readers ignore this information and might give wrong results.

Reading simple linear 1D WCS
----------------------------

One of the most common and simple ways that dispersion is encoded in FITS files is linear dispersion using four keywords::

    CRVAL1  =    4402.538203477947 # Dispersion at reference pixel
    CRPIX1  =                    1 # Reference pixel
    CDELT1  =      1.3060348033905 # Dispersion by pixel
    CUNIT1  = 'Angstrom'


One can easily create a simple wcs file from this information::

    >>> from specutils.wcs import specwcs
    >>> from astropy.io import fits
    >>> from astropy import units as u
    >>> header = fits.getheader('myfile.fits')
    >>> dispersion_start = header['CRVAL1'] - (header['CRPIX1'] - 1) * header['CDELT1']
    >>> linear_wcs = specwcs.Spectrum1DPolynomialWCS(degree=1, c0=dispersion_start,
    >>>                                              c1=header['CDELT1'],
    >>>                                              unit=u.Unit(header['CUNIT1']))
    >>> flux = fits.getdata('myfile.fits')
    >>> myspec = Spectrum1D(flux=flux, wcs=linear_wcs)

However, this process requires too much code for the user. As this is such a common format, there exists a reader that generates a
Spectrum1D object from a fits file directly ::

   >>> from specutils.io import read_fits
   >>> myspec = read_fits.read_fits_spectrum1d('myfile.fits')

`read_fits_spectrum1d` is the only method required to read various types of FITS formats listed in `this page <http://iraf.net/irafdocs/specwcs.php>`_.
This method detects the type of WCS, and uses the correct reader accordingly.

Currently only linear one-dimensional WCS is implemented, but the examples should give a guide to implement more complex
or WCS created from different keywords.


Reading FITS "multispec" format WCS
------------------------------------

The multispec format holds multiple one-dimensional spectra in a single file. The current 1-D readers can be used to read
such a file. A list of spectra is returned whenever the input file is a multispec file.
Here is an example of reading a simple FITS multispec format::

    >>> from specutils.io import read_fits
    >>> spectra_list = read_fits.read_fits_spectrum1d('mymultispec.fits')
    >>> spectra_list[0]
    Spectrum1D([ 338.06109619,  395.59234619,  326.0012207 , ...,
                 660.0098877 ,  686.54498291,  689.58374023])
    >>> spectra_list[1].dispersion
    <Quantity [ 8293.44875263, 8293.40459999, 8293.36044556,...,
                8166.53073537, 8166.48250242, 8166.43426803] Angstrom>

The multispec format supports various functions to map the pixel indices to dispersion. Currently the following multispec formats
are supported:

 * Linear Dispersion Function
 * Log-linear Dispersion Function
 * Legendre Polynomial Dispersion Function
 * Chebyshev Polynomial Dispersion Function
 * Linear Spline Dispersion Function
 * Cubic Spline Polynomial Dispersion Function

.. warning::
    Linear and cubic spline functions implementations haven't been tested due to lack of test
    files. If you have one of these files, please post at `astropy-dev <https://groups.google.com/forum/#!forum/astropy-dev>`_.


Reference/API
-------------

.. automodapi:: specutils.io.read_fits
    :no-inheritance-diagram:
    :skip: OrderedDict, Spectrum1D
