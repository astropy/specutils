===================
Spectrum Fitting
===================

Specutils has the ability to fit the flux of a `~specutils.Spectrum1D` object.
The fit routine takes the `~specutils.Spectrum1D` object and a list of 
`~astropy.modeling.Model` that have initial guesses for each of the parameters. 

The internal functionality uses `~astropy.modeling.fitting` routines.  The flux
is extracted from the `~specutils.Spectrum1D` and is passed along with a 
compound model created from the model initial guesses.

Model Fitting
-------------

The first step is to create a set of models with initial guesses as the parameters. Even
better is to include a set of bounds for each parameter, but that is optional. 

.. code-block:: python

    >>> from astropy.modeling import models
    >>> from specutils import Spectrum1D
    >>> from specutils import fitting
    >>> import astropy.units as u
    >>> import numpy as np
    >>> np.random.seed(42)


For the purpose of the example, build a ``spectrum1d`` variable that will be used in the fitting:

.. code-block:: python

	>>> # Create the wavelength array
    >>> wave_um = np.linspace(0.4, 1.05, 100)

	>>> # Create the models
    >>> g1 = models.Gaussian1D(amplitude=2000, mean=0.52, stddev=0.01)
    >>> g2 = models.Gaussian1D(amplitude=500, mean=0.64, stddev=0.02)
    >>> g3 = models.Gaussian1D(amplitude=-350, mean=0.78, stddev=0.01)
    >>> ramp = models.Linear1D(slope=300, intercept=0.0)

	>>> # Create the flux array
    >>> base_flux = g1(wave_um) + g2(wave_um) + g3(wave_um) + ramp(wave_um)

	>>> # Create the `specutils.Spectrum1D` object.
    >>> flux_e1 = base_flux + 100*np.random.random(base_flux.shape)
    >>> spectrum1d = Spectrum1D(spectral_axis=wave_um*u.um, flux=flux_e1*u.mJy)


.. figure:: img/fitting_original.jpg
   :figwidth: 500px
   :align: center

   Spectrum to be fit.

Now that there is a ``spectrum1d`` to fit, the real fitting setup must happen.  First create
the `astropy.modeling.Model` to be used in the fitting routine.

.. code-block:: python

    >>> # Create the initial guesses
    >>> fg1 = models.Gaussian1D(amplitude=1900, mean=5175, stddev=200,
            bounds = {'amplitude': (1700, 2400), 
                      'mean': (4800, 5500), 
                      'stddev': (1.1754943508222875e-38, 500)})
    >>> fg2 = models.Gaussian1D(amplitude=400, mean=6232, stddev=150,
            bounds = {'amplitude': (200, 500), 
                      'mean': (5500, 7000), 
                      'stddev': (1.1754943508222875e-38, 300)})
    >> fg3 = models.Gaussian1D(amplitude=-500, mean=7900, stddev=150,
            bounds = {'amplitude': (-800, -300), 
                      'mean': (7300, 8200), 
                      'stddev': (1.1754943508222875e-38, 300)})
    >>> framp = models.Linear1D(slope=320/10000, intercept=0)

Now comes the actual fitting:

.. code-block:: python

    >>> # And... now we can do the fitting
    >>> fitted_models = fitting.fitmodels(spectrum1d, [fg1, fg2, fg3, framp])

The output of the fitting routine is another set of models whose parameters contain the result of the fit.
(There is a corresponding output model for each input model.)

.. figure:: img/fitting_original_fit.jpg
   :figwidth: 500px
   :align: center

   Original spectrum (blue) and fitted spectrum (orange).

.. automodapi:: specutils.fitting
