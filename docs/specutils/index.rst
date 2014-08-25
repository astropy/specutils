Spectrum Utilities
==================

The ``specutils``-package is intended to implement utilities and data structures to handle astronomical spectra.
The initial goal is to focus on electromagnetic spectra, but the end goal is to support spectra of any kind.

Most commonly, astronomical spectra are given in an array with one to three dimensions.
The `~astropy.nddata.NDData`-class is an ideal storage format for this data, as
it provides the possibility to store the array but also various forms of meta-data.

Transformation functions that map the pixel indices to physical meaningful coordinates
(in a one-dimensional case - e.g. Angstrom) are stored as metadata in a `wcs`-keyword (WCS meaning World Coordinate System transform).
More information about Spectral WCS can be found :doc:`here <specwcs>`. The `flag`-keyword is designed to hold pixel specific
meta information (i.e. sky arrays). The `mask`-keyword is a special kind of flag that only carries boolean values and is used to
create `~np.ma.MaskedArrays` out of NDData objects. Finally, the `meta`-keyword is designed to hold general metadata in a dictionary-like
format. These metadata options comprise the most needed kinds of meta data.

The package is subdivided into multiple parts:

Definition of 1-, 2- and 3-D spectral containers including transformations from pixel coordinates to
physical coordinates (World Coordinate System - WCS).

**Spectral Containers**

* :doc:`Spectrum1D <spectrum1d>`
* :doc:`Spectrum2D <spectrum2d>`
* :doc:`Spectrum3D <spectrum3d>`

.. note::
    Currently, only 1-dimensionanal spectra are supported

**World Coordinate System**

* :doc:`Spectral World Coordinate Systems <specwcs>`

**I/O**

* :doc:`Reading spectra from FITS files <read_fits>`
* :doc:`Writing spectra to FITS files <write_fits>`

**Spectral Tools and Utilities**

* :doc:`Interstellar Extinction <extinction>`



