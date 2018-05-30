===================
Spectrum Fitting
===================

Specutils has the ability to fit the flux of a `specutils.spectra.Spectrum1D` object.
The fit routine takes the `specutils.spectra.Spectrum1D` object and a list of 
`astropy.modeling.models` that have initial guesses for each of the parameters. 

The internal functionality uses `astropy.modeling.fitting` routines.  The flux
is extracted from the `specutils.spectra.Spectrum1D` and is passed along with a 
compound model created from the model initial guesses.

Basic Fitting
-------------

The first step is to create a set of models with initial guesses as the parameters. Even
better is to include a set of bounds for each parameter, but that is optional. 

.. code-block:: python

    >>> from astropy.modeling import models
    >>> from specutils import Spectrum1D
    >>> from specutils import fitting
    >>> import astropy.units as u
    >>> import numpy as np


For the purpose of the example, build a spectrum1d variable that will be used in the fitting::

	>>> # Create the wavelength array
    >>> wave_um = np.linspace(0.4, 1.05, 100)

	>>> # Create the models
    >>> g1 = models.Gaussian1D(amplitude=2000, mean=0.56, stddev=0.01)
    >>> g2 = models.Gaussian1D(amplitude=500, mean=0.62, stddev=0.02)
    >>> g3 = models.Gaussian1D(amplitude=-350, mean=0.52, stddev=0.01)
    >>> ramp = models.Linear1D(slope=300, intercept=0.0)

	>>> # Create the flux array
    >>> base_flux = g1(wave_um) + g2(wave_um) + g3(wave_um) + ramp(wave_um)

	>>> # Create the `specutils.spectra.Spectrum1D` object.
    >>> flux_e1 = base_flux + 200*np.random.random(base_flux.shape)
    >>> specrum1d = Spectrum1D(spectral_axis=wave_um*u.um, flux=flux_e1*u.mJ)

Now that there is a spectrum1d to fit, the real fitting setup must happen.  First create
the `astropy.modeling.models` to be used in the fitting routine.

    >>> # Create the initial guesses
    >>> fg1 = models.Gaussian1D(amplitude=1900, mean=5575, stddev=0.02*1000,
            bounds = {'amplitude': (1700, 2100), 'mean': (5000, 6000), 'stddev': (1.1754943508222875e-38, 200)})
    >>> fg2 = models.Gaussian1D(amplitude=400, mean=6232, stddev=0.01*1000,
            bounds = {'amplitude': (200, 600), 'mean': (5500, 7000), 'stddev': (1.1754943508222875e-38, 200)})
    >> fg3 = models.Gaussian1D(amplitude=-300, mean=8005, stddev=0.01*1000,
            bounds = {'amplitude': (-500, -100), 'mean': (7500, 8500), 'stddev': (1.1754943508222875e-38, 200)})
    >>> framp = models.Linear1D(slope=0.8, intercept=0)

Now comes the actual fitting::

    >>> # And... now we can do the fitting
    >>> fitted_models = fitting.fit_models(simulated_spectra.s1_um_mJy_e1, [fg1, fg2, fg3, framp])

The output of the fitting routine is another set of models whose parameters contain the result of the fit.
(There is a corresponding output model for each input model.)
