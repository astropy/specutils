.. _doc_installation:

Installation
============

Dependencies
------------

In the future, SpecViz will be distributed through package managers like
`Anaconda <https://anaconda.org>`_ and `Homebrew <http://brew.sh>`_ that will obviate the need for manual installation from
source. Most of these will be handled automatically by the setup functions,
with the exception of ``PyQt``/``PySide``.

* Python 3 (recommended) or Python 2
* PyQt5 (recommended), PyQt4, or PySide
* Astropy
* Numpy
* Scipy
* PyQtGraph


Installing PyQt/PySide
----------------------
The easiest way to install PyQt/PySide is through some package manager.
Please keep in mind that PyQt5 is the recommended PyQt implementation as
`Qt4 development and support has ended <http://blog.qt.io/blog/2015/05/26/qt-4-8-7-released/>`_.

Below are instructions for installing using *either* Homebrew *or* Anaconda.

PyQt5
^^^^^

Homebrew
   `Install using Homebrew for Qt5 <http://brewformulas.org/Pyqt5>`_.

Anaconda
   Installing PyQt5 with Anaconda will require installing from the Spyder-IDE
   channel as it is not currently a core package (but they're working on it).

   Further, Anaconda will panic if you have both PyQt4 and PyQt5 installed in
   the same environment. To work around this, it is strongly suggested you
   simply create a new virtual environment and install PyQt5 there::

    $ conda create -n pyf_env python=3.5
    $ source activate pyf_env
    $ conda install --channel https://conda.anaconda.org/spyder-ide pyqt5

.. note::
   PyQt5 **does not** require Python 3. If you wish, you can create your
   virtual environment using Python 2 by specifying the version as shown above
   (e.g. ``python=2.7``).

PyQt4
^^^^^

Homebrew
   `Install using Homebrew for Qt4 <http://brewformulas.org/Pyqt4>`_.

Anaconda
   Install using Anaconda::

    $ conda install pyqt


Installing SpecViz
------------------

Clone the SpecViz repository somewhere on your system, and install locally using
``pip``. If you are using an Anaconda virtual environment, please be sure to
activate it first before installing: ``$ source activate pyf_env``.

::

    $ git clone https://github.com/nmearl/specviz.git
    $ cd specviz
    $ git checkout tags/v0.1rc2
    $ pip install -r requirements.txt

.. note::

   This uses the ``pip`` installation system, so please note that

   1. You need to have ``pip`` installed (included in most Python
      installations).
   2. You do **not** need to run ``python setup.py install``.
   3. You do **not** need to install the dependencies by hand (except for PyQt).

   Likewise, the ``pip`` command will use your default Python to install.
   You can specify by using ``pip2`` or ``pip3``, if you're not using a virtual
   environment.


Known Issues
------------

On a Mac with Qt5, depending on exactly how you have set up Anaconda, you might
see the following error after following the above instructions::

    This application failed to start because it could not find or load the Qt platform plugin "cocoa".

    Reinstalling the application may fix this problem.

If you see this message, you have encountered an incompatibility between
Anaconda's packaging of Qt4 and Qt5. The workaround is to uninstall Qt4 with the
following command::

    $ conda uninstall pyqt qt

and SpecViz should now happily run.


.. _doc_launching:

Launching SpecViz
=================

Once PyQt/PySide and the other dependencies are installed, SpecViz can be
launched from the command line::

    $ specviz

Again, if you're using an Anaconda virtual environment, please be sure to
activate it before launching SpecViz.

If you only wish to inspect a single FITS or ASCII file using the default
:ref:`doc_custom_loaders` file formatting, you can also pass in the filename
as a command line argument, as follows::

    $ specviz filename
