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
    systems is in development. Currently, FITS readers ignore this information and might give wrong results.

Reading simple linear 1D WCS
----------------------------

One of the most common and simple ways that dispersion is encoded in FITS files is linear dispersion using four keywords::

    CRVAL1  =    4402.538203477947
    CRPIX1  =                    1
    CDELT1  =      1.3060348033905
    CUNIT1  = 'Angstrom'

One can easily create a simple wcs file from this information::

    >>> from specutils.wcs import specwcs
    >>> from astropy.io import fits
    >>> from astropy import units as u
    >>> header = fits.getheader('myfile.fits')
    >>> dispersion_start = dispersion_start - (header['CRPIX1'] - 1)  * dispersion_delta
    >>> linear_wcs = specwcs.Spectrum1DPolynomialWCS(degree=1, c0=dispersion_start,
                                                        c1=header['CDELT1'], unit=u.Unit(header['CUNIT1']))
    >>> flux = fits.getdata('myfile.fits')
    >>> myspec = Spectrum1D(flux=flux, wcs=linear_wcs)

As this is such a common format, there is a WCS reader available which generates a linear_wcs::

    >>> from astropy.io import fits
    >>> from specutils.io import read_fits
    >>> fits_wcs_info = read_fits.FITSWCSSpectrum(fits.getheader('myfile.fits')) # parsing the FITS WCS information
    >>> linear_wcs = read_fits_wcs_linear1d(fits_wcs_info)

Finally, there exists a reader that generates a Spectrum1D object from a fits file. The reader will iterate over the
available FITS 1D WCS readers until a matching one is found::

   >>> from specutils.io import read_fits
   >>> myspec = read_fits.read_fits_spectrum1d('myfile.fits'):

Currently only linear one-dimensional WCS is implemented, but the examples should give a guide to implement more complex
or WCS created from different keywords.


Reading FITS "multispec" format WCS
------------------------------------

Here is an example of reading a simple FITS multispec format. The output will not be a two-dimensional Spectrum, but
a list of `~specutils.Spectrum1D` objects::

    >>> from specutils.io import read_fits
    >>> read_fits.read_fits_multispec_to_list('mymultispec.fits')

Internally, the function again uses the `~FITSWCSSpectrum` parser object.



Reference/API
-------------

.. automodapi:: specutils.io.read_fits
    :no-inheritance-diagram:
    :skip: OrderedDict, Spectrum1D
