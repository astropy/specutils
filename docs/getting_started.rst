****************
Creating Spectra
****************

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


Spectra with units
------------------

It's also possible to include units. This can be done either by passing in
:class:`astropy.units.Quantity` arrays, or by specifying the units explicitly.

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt
    import astropy.units as u

    from specutils import Spectrum1D

    flux = np.random.sample(200)
    wave = np.arange(1100, 1300)

    # Specifying units explicitly
    spec1d = Spectrum1D(flux, spectral_axis=wave, unit=u.Jy, spectral_axis_unit=u.AA)

    # Using astropy quantities
    spec1d = Spectrum1D(flux * u.Jy, spectral_axis=wave * u.AA)

    plt.plot(spec1d.wavelength, spec1d.flux)

.. image:: img/quick_start2.png


Defining WCS
------------