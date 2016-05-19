.. SpecViz documentation master file, created by
   sphinx-quickstart on Mon Feb  8 02:58:02 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SpecViz: 1D Spectral Visualization Tool
=======================================

SpecViz is a tool for visualization and quick-look analysis of 1D astronomical
spectra. It is written in the Python programming language, and therefore can be
run anywhere Python is supported (see :ref:`doc_installation`). SpecViz is
capable of reading data from FITS and ASCII tables (see :ref:`doc_custom_loaders`).

SpecViz allows spectra to be easily plotted and examined. It supports
instrument-specific data quality handling, flexible spectral units conversions,
custom plotting attributes, plot annotations, tiled plots, and other features.

SpecViz notably includes a measurement tool for spectral lines which
enables the user, with a few mouse actions, to perform and record measurements.
It has a model fitting capability that enables the user to create simple
(e.g., single Gaussian) or multi-component models (e.g., multiple Gaussians for
emission and absorption lines in addition to regions of flat continuua).
SpecViz incorporates various methods for fitting such models to data. For more
details, see :ref:`doc_model_fitting`.

Furthermore, SpecViz allows for overplotting or simple combining of spectra at
different wavelengths.

SpecViz will soon include the ability to
   - Process spectra using averaging, splicing, detrending, and Fourier filtering tools.
   - Support overplotting and interactively renormalizing data from spectral templates.
   - Overplot of spectral line identifications taken from a variety of line lists.
   - And more...

Demo
----

.. image:: https://i.vimeocdn.com/video/571749719_640.jpg
   :target: https://vimeo.com/167441251


Installation and Setup
----------------------

.. toctree::
   :maxdepth: 2

   installation


Using SpecViz
-------------

.. toctree::
   :maxdepth: 2

   viewer
   model_fitting
   custom_loaders


References/API
--------------

.. toctree::
   :maxdepth: 1

   api/index


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
