Reading Spectra from a FITS file
=================================

FITS is one of the most common ways to store spectra. In most cases the flux values are stored as the pixel values
in a 1-, 2- or 3-D arrays. The mapping of axes indices to physically meaningful values (e.g. to wavelength) is encoded
in so FITS keywords. This document describes how to extract this information and read it in to the various spectral
classes. A lot of information about the storage is taken from
`The IRAF/NOAO Spectral World Coordinate Systems <http://iraf.net/irafdocs/specwcs.php>`_.

Logical/physical

Reading simple linear 1D WCS
----------------------------

One of the most common and simple
