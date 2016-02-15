.. Pyfocal documentation master file, created by
   sphinx-quickstart on Mon Feb  8 02:58:02 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Pyfocal: 1D Spectral Visualization Tool
=======================================

Pyfocal is a tool for 1-D spectral visualization and analysis of astronomical
spectrograms. It is written in Python thus can be run anywhere Python is
supported (see :ref:`doc_installation`).
Pyfocal is capable of reading data from FITS and ASCII tables
(see :ref:`doc_custom_loaders`).

Once ingested, data can be plotted and examined with a large selection of
custom settings. Pyfocal supports instrument-specific data quality handling,
flexible spectral units conversions, custom plotting attributes,
plot annotations, tiled plots, etc.

A spectral feature quick measurement tool enables the user, with a few mouse
actions, to perform and record measurements on selected spectral features.

Pyfocal can be used to build wide-band SEDs, overploting or combining data from
the same astronomical source taken with different instruments and/or spectral
bands. Data can be further processed with averaging, splicing, detrending,
and Fourier filtering tools.

Pyfocal has a spectral model fitting capability that enables the user to work
with multi-component models in a number of ways, and fit models to data.
For more details, see :ref:`doc_model_fitting`.

Support exists for overplotting and interactively renormalize data from
spectral templates.

Pyfocal can overplot spectral line identifications taken from a variety of
line lists.


.. _doc_installation:

Installation
------------

.. toctree::
   :maxdepth: 2

   installation


Using Pyfocal
-------------

.. toctree::
   :maxdepth: 2

   model_fitting
   custom_loaders
   api/index


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
