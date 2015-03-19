specutils
=========

This package contains an Astropy affiliated package for operations and tools on astronomical spectra.
It is subdivided into multiple parts:

Spectral Containers
-------------------------
Spectrum1D | Spectrum2D | Spectrum3D

    **Note:**
    Currently, only 1-dimensional spectra are currently supported

I/O
-------------------------
* Reading FITS *simple linear 1D WCS*
**read_fits_spectrum1d** is the only method required to read various types of FITS formats listed in `this page`_. This method detects the type of WCS, and uses the correct reader accordingly.

.. code-block:: python

    from specutils.io import read_fits
    myspec = read_fits.read_fits_spectrum1d('myfile.fits')

Currently only linear one-dimensional WCS is implemented, but the examples should give a guide to implement more complex or WCS created from different keywords.

* Reading FITS *multispec format WCS*
The multispec format holds multiple one-dimensional spectra in a single file. The current 1-D readers can be used to read such a file. A list of spectra is returned whenever the input file is a multispec file.

.. code-block:: python

    from specutils.io import read_fits
    spectra_list = read_fits.read_fits_spectrum1d('mymultispec.fits')
    spectra_list[0]
    Spectrum1D([ 338.06109619,  395.59234619,  326.0012207 , ...,
                 660.0098877 ,  686.54498291,  689.58374023])
    spectra_list[1].dispersion
    <Quantity [ 8293.44875263, 8293.40459999, 8293.36044556,...,
                8166.53073537, 8166.48250242, 8166.43426803] Angstrom>

* Writing spectra to FITS files
**write_fits.write** method deciphers the type of object passed automatically and writes spectra to the given file in FITS format. 

.. code-block:: python

    from specutils.io import read_fits, write_fits
    myspec = read_fits.read_fits_spectrum1d('myfile.fits')
    write_fits.write(myspec, 'mynewfile.fits')

Writing to a file in multispec FITS format is identical.

World Coordinate System
-------------------------
* Spectral 1D WCS
The simplest WCS for 1D is a lookup table. This is often found in a ascii tables where one column is the lookup table (wavelength array) and one column is the flux array. In terms of the functional transform view: the lookup table represents a parameter (c0):

.. code-block:: python

    from specutils import Spectrum1D
    from specutils.wcs import specwcs
    wave = numpy.arange(6000, 9000) # wavelength
    flux = numpy.random.random(3000) # flux
    lookup_table_wcs = specwcs.Spectrum1DLookupWCS(wave, unit='angstrom') # Creating the wcs
    lookup_table_wcs(0) # doing the transformation for pixel 0
    <Quantity 6000.0 Angstrom>
    Spectrum1D(flux=flux, wcs=lookup_table_wcs)
    Spectrum1D([ 0.66118716,  0.39584688,  0.81210479, ...,  0.5238203 ,
                 0.05986459,  0.11621466])

For more information about the Spectrum1d object go to `Spectrum 1D`_.

Another common WCS is the **linear dispersion** and commonly serialized (encoded) to FITS keyword headers. For linear dispersion we are using the general `Spectrum1DPolynomialWCS`_ WCS.

The `CompositeWCS`_ and `WeightedCombinationWCS`_ models can be useful to combine different WCS.
Another important model available is the `DopplerShift`_ model. This model is specifically for calculating the doppler shift from velocity (v).

In addition, the following WCS models exist as well:
    * `Spectrum1DIRAFLegendreWCS`_
    * `Spectrum1DIRAFChebyshevWCS`_
    * `Spectrum1DIRAFBSplineWCS`_
    * `MultispecIRAFCompositeWCS`_
These models are just IRAF extensions of already existing models. This extensions enable the user to read and write from IRAF FITS files.

Spectral Tools and Utilities
--------------------------------
* Interstellar Extinction calculations
This module contains extinction law functions. See the documentation for individual functions.

`Full Documentation`_ 

.. image:: https://travis-ci.org/astropy/specutils.png?branch=master
  :target: https://travis-ci.org/astropy/specutils

.. image:: https://coveralls.io/repos/astropy/specutils/badge.png
  :target: https://coveralls.io/r/astropy/specutils

.. _Full Documentation: http://astroquery.readthedocs.org
.. _this page: http://iraf.net/irafdocs/specwcs.php
.. _Spectrum 1D: http://specutils.readthedocs.org/en/latest/specutils/spectrum1d.html
.. _Spectrum1DPolynomialWCS: http://specutils.readthedocs.org/en/latest/api/specutils.wcs.specwcs.Spectrum1DPolynomialWCS.html#specutils.wcs.specwcs.Spectrum1DPolynomialWCS
.. _CompositeWCS: http://specutils.readthedocs.org/en/latest/api/specutils.wcs.specwcs.CompositeWCS.html#specutils.wcs.specwcs.CompositeWCS
.. _WeightedCombinationWCS: http://specutils.readthedocs.org/en/latest/api/specutils.wcs.specwcs.WeightedCombinationWCS.html#specutils.wcs.specwcs.WeightedCombinationWCS
.. _DopplerShift: http://specutils.readthedocs.org/en/latest/api/specutils.wcs.specwcs.DopplerShift.html#specutils.wcs.specwcs.DopplerShift
.. _Spectrum1DIRAFLegendreWCS: http://specutils.readthedocs.org/en/latest/api/specutils.wcs.specwcs.Spectrum1DIRAFLegendreWCS.html#specutils.wcs.specwcs.Spectrum1DIRAFLegendreWCS
.. _Spectrum1DIRAFChebyshevWCS: http://specutils.readthedocs.org/en/latest/api/specutils.wcs.specwcs.Spectrum1DIRAFChebyshevWCS.html#specutils.wcs.specwcs.Spectrum1DIRAFChebyshevWCS
.. _Spectrum1DIRAFBSplineWCS: http://specutils.readthedocs.org/en/latest/api/specutils.wcs.specwcs.Spectrum1DIRAFBSplineWCS.html#specutils.wcs.specwcs.Spectrum1DIRAFBSplineWCS
.. _MultispecIRAFCompositeWCS: http://specutils.readthedocs.org/en/latest/api/specutils.wcs.specwcs.MultispecIRAFCompositeWCS.html#specutils.wcs.specwcs.MultispecIRAFCompositeWCS
