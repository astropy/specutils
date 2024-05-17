Specutils
=========

|CI| |Coverage| |Docs| |Astropy|

.. |CI| image:: https://github.com/astropy/specutils/workflows/CI/badge.svg
   :target: https://github.com/astropy/specutils/actions
   :alt: GitHub Actions CI Status

.. |Coverage| image:: https://codecov.io/github/astropy/specutils/branch/main/graph/badge.svg
   :target: https://codecov.io/github/astropy/specutils
   :alt: Coverage Status

.. |Docs| image:: https://readthedocs.org/projects/specutils/badge/?version=latest
   :target: https://specutils.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. |Astropy| image:: https://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
   :target: https://www.astropy.org/
   :alt: Powered by Astropy

Specutils is an `Astropy affiliated package <http://affiliated.astropy.org/>`_
with the goal of providing a shared set of Python representations of
astronomical spectra and basic tools to operate on these spectra. The effort is
also meant to be a "hub", helping to unite the Python astronomical spectroscopy
community around shared effort, much as Astropy is meant to for the wider
astronomy Python ecosystem. This broader effort is outlined in the
`APE13 document <https://github.com/astropy/astropy-APEs/blob/main/APE13.rst>`_.

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

This project is Copyright (c) Specutils Developers and licensed under
the terms of the BSD 3-Clause license. This package is based upon
the `Astropy package template <https://github.com/astropy/package-template>`_
which is licensed under the BSD 3-clause license. See the ``licenses`` folder for
more information.

Contributing
------------

We love contributions! specutils is open source,
built on open source, and we'd love to have you hang out in our community.

**Imposter syndrome disclaimer**: We want your help. No, really.

There may be a little voice inside your head that is telling you that you're not
ready to be an open source contributor; that your skills aren't nearly good
enough to contribute. What could you possibly offer a project like this one?

We assure you - the little voice in your head is wrong. If you can write code at
all, you can contribute code to open source. Contributing to open source
projects is a fantastic way to advance one's coding skills. Writing perfect code
isn't the measure of a good developer (that would disqualify all of us!); it's
trying to create something, making mistakes, and learning from those
mistakes. That's how we all improve, and we are happy to help others learn.

Being an open source contributor doesn't just mean writing code, either. You can
help out by writing documentation, tests, or even giving feedback about the
project (and yes - that includes giving feedback about the contribution
process). Some of these contributions may be the most valuable to the project as
a whole, because you're coming to the project with fresh eyes, so you can see
the errors and assumptions that seasoned contributors have glossed over.

Note: This disclaimer was originally written by
`Adrienne Lowe <https://github.com/adriennefriend>`_ for a
`PyCon talk <https://www.youtube.com/watch?v=6Uj746j9Heo>`_, and was adapted by
specutils based on its use in the README file for the
`MetPy project <https://github.com/Unidata/MetPy>`_.

If you locally cloned this repo before 22 Mar 2021
--------------------------------------------------

The primary branch for this repo has been transitioned from ``master`` to ``main``.  If you have a local clone of this repository and want to keep your local branch in sync with this repo, you'll need to do the following in your local clone from your terminal::

   git fetch --all --prune
   # you can stop here if you don't use your local "master"/"main" branch
   git branch -m master main
   git branch -u origin/main main

If you are using a GUI to manage your repos you'll have to find the equivalent commands as it's different for different programs. Alternatively, you can just delete your local clone and re-clone!
