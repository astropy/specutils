*********
Specutils
*********

The specutils package provides a basic interface for the loading, manipulation,
and low-level analysis of spectroscopic data. The intention is for the
generic data containers and accompanying modules to provide a basis for the
astronomical community upon which more robust tooling can be built.

The basic data container of specutils is the :class:`~specutils.spectra.Spectrum1D`
class. It contains logic to handle multi-dimensional flux data, spectral axis
in various forms (wavelenth, frequency, energy, velocity, etc.), convenient and
unobtrusive wcs support, and uncertainty handling. This core container is
supported by a few other features including unit support, custom data loading,
arithmetic operation support, and multi-dimensional flux data support

.. note:: support for collections of spectra with distinct spectral axes
          information is slotted for a future release.

Quick Start
-----------

Generating a spectrum object is quite easy. The simplest way would look
something like

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt

    from specutils import Spectrum1D

    flux = np.random.sample(200)
    wave = np.arange(1100, 1300)

    spec1d = Spectrum1D(flux, spectral_axis=wave)

    plt.plot(spec1d.wavelength, spec1d.flux)

.. image:: img/quick_start.png