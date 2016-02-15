.. _doc_installation:

Installation
============

Dependencies
------------
In the future, Pyfocal will be distributed through package managers like Anaconda and Homebrew that will obviate the
need for manual installation from source. Most of these will be handled automatically by the setup functions,
with the exception of ``PyQt``/``PySide``.

* Python 3 (recommended), or Python 2
* PyQt5 (recommended), PyQt4, or PySide
* Astropy
* Numpy
* Scipy
* PyQtGraph

Installing PyQt/PySide
^^^^^^^^^^^^^^^^^^^^^^
The easiest way to install PyQt/PySide is through some package manager. Please keep in
mind that PyQt5 is the
recommended PyQt implementation as `Qt4 development and support has ended <http://blog.qt
.io/blog/2015/05/26/qt-4-8-7-released/>`_.

Below are instructions for installing using *either* Homebrew *or* Anaconda.

For PyQt5
"""""""""

Homebrew:
   `Install using Homebrew for Qt5 <http://brewformulas.org/Pyqt5>`_.

Anaconda:
   Installing PyQt5 with Anaconda will require installing from the Spyder-IDE channel as it is not currently a core
   package (but they're working on it).

   Further, Anaconda will panic if you have both PyQt4 and PyQt5 installed in the same environment. To work around
   this, it is strongly suggested you simply create a new virtual environment and install PyQt5 there::

    $ conda create -n pyf_env python=3.5
    $ source activate pyf_env
    $ conda install --channel https://conda.anaconda.org/spyder-ide pyqt5

.. note::
   PyQt5 **does not** require Python 3. If you wish, you can create your virtual environment using Python 2 by
   specifying the version as shown above (e.g. ``python=2.7``).

For PyQt4
"""""""""

Homebrew
   `Install using Homebrew for Qt4 <http://brewformulas.org/Pyqt4>`_.

Anaconda
   Install using Anaconda: ``conda install pyqt``.

Installing from source
----------------------
Go ahead and clone the repository somewhere on your system, and install locally using ``pip``. If you are using an
Anaconda virtual environment, please be sure to activate it first before installing: ``$ source
activate pyf_env``.

::

    $ git clone https://github.com/nmearl/pyfocal.git
    $ cd pyfocal
    $ git checkout tags/v0.1rc2
    $ pip install -r requirements.txt

.. note::

   This uses the ``pip`` installation system, so please note that

   1. you need to have ``pip`` installed (included in most Python installations),
   2. you do **not** need to run ``python setup.py install``,
   3. you do **not** need to install the dependencies by hand (except for PyQt).

   Likewise, the ``pip`` command will use your default python to install. You can specify by using ``pip2`` or ``pip3``, if you're not using e.g. a virtual environment.

Launching from command line
---------------------------
Once PyQt/PySide and the other dependencies are installed, Pyfocal can be launched from the command line ::

    $ pyfocal

Again, if you're using an Anaconda virtual environment, please be sure to activate before launching Pyfocal.
