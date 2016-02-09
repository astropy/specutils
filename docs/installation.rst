.. _`Installation`:

Installation
============

Dependencies
------------
In the future, Pyfocal will be distributed through package managers like Anaconda and Homebrew that will obviate the
need for manual installation of the dependencies. Most of these will be handled automatically by the setup functions,
with the exception of ``PyQt``/``PySide``.

* Python 3 (recommended), or Python 2
* PyQt5 (recommended), PyQt4, or PySide
* Astropy
* Numpy
* PyQtGraph
* py_expression_eval
* qtpy

Installing PyQt/PySide
^^^^^^^^^^^^^^^^^^^^^^
The easiest way to install PyQt/PySide is through some package manager. Please keep in mind that PyQt5 is the
recommended PyQt implementation as `Qt4 development and support has ended <http://blog.qt
.io/blog/2015/05/26/qt-4-8-7-released/>`_.

**For PyQt5**

* `Install using Homebrew <http://brewformulas.org/Pyqt5>`_.
* PyQt5 is not a core package yet for Anaconda (but they promise by the beginning of 2016, it will be). Currently, to install with Anaconda, you should tap the Spyder-IDE channel::

    $ conda install --channel https://conda.anaconda.org/spyder-ide pyqt5

**For PyQt4**

* `Install using Homebrew <http://brewformulas.org/Pyqt4>`_
* Install using Anaconda: ``conda install pyqt``

Installing from source
----------------------
Go ahead and clone the repository somewhere on your system, and install locally using ``pip``.

.. note::

   This uses the ``pip`` installation system.

   1. You need to have ``pip`` installed
   2. You do **not** need to run ``python setup.py install``

   Likewise, the ``pip`` command will use your default python to install. You can specify by using ``pip2`` or ``pip3``, if you're not using e.g. a virtual environment.

::

    $ git clone https://github.com/nmearl/pyfocal.git
    $ cd pyfocal
    $ pip install -r requirements.txt

Launching from command line
---------------------------
Once PyQt/PySide and the other dependencies are installed, Pyfocal can be launched from the command line ::

    $ pyfocal

