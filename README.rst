Specutils
=========

.. image:: https://travis-ci.org/astropy/specutils.svg?branch=master
    :target: https://travis-ci.org/astropy/specutils

.. image:: https://readthedocs.org/projects/specutils/badge/?version=latest
    :target: http://specutils.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org/

Specutils is an `Astropy affiliated package <http://affiliated.astropy.org/>`_
with the goal of providing a shared set of Python representations of
astronomical spectra and basic tools to operate on these spectra. The effort is
also meant to be a "hub", helping to unite the Python astronomical spectroscopy
community around shared effort, much as Astropy is meant to for the wider
astronomy Python ecosystem. This broader effort is outlined in the
`APE13 document <https://github.com/astropy/astropy-APEs/blob/master/APE13.rst>`_.

Note that Specutils is not intended to meet all possible spectroscopic analysis or
reduction needs. While it provides some standard analysis functionality
(following the  Python philosophy of "batteries included"), it is best thought
of as a "tool box" that provides pieces of functionality that can be used to
build a particular scientific workflow or higher-level analysis tool.  To that
end, it is also meant to facilitate connecting together disparate reduction
pipelines and analysis tools through shared Python representations of
spectroscopic data.

Getting Started with Specutils
------------------------------

For details on installing and using Specutils, see the
`specutils documentation <http://specutils.readthedocs.io/en/latest/>`_.

Note that Specutils now only supports Python 3. While some functionality may
continue to work on Python 2, it is not tested and support cannot be guaranteed
(due to the sunsetting of Python 2 support by the Python and Astropy development
teams).

License
-------

Specutils is licensed under a 3-clause BSD style license. Please see the LICENSE.rst file.
