****************
Creating Spectra
****************

Spectra With Units
------------------

It's also possible to include units.

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt
    import astropy.units as u

    from specutils import Spectrum1D

    flux = np.random.sample(200) * u.Jy
    wave = np.arange(1100, 1300) * u.AA

    spec1d = Spectrum1D(flux, spectral_axis=wave)

    plt.plot(spec1d.wavelength, spec1d.flux)

.. image:: img/quick_start2.png


Defining WCS
------------