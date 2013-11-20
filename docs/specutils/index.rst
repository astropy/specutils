Spectrum Utilities
==================

The ``specutils``-package is intended to implement utilities and data structures to handle astronomical spectra.
The initial goal is to focus on electromagnetic spectra, but the end goal is to support spectra of any kind.

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



**Spectral Tools and Utilities**

* :doc:`Interstellar Extinction <extinction>`



