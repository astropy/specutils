.. _doc_model_fitting:

Model Fitting
=============

Pyfocal utilizes
`Astropy Models and Fitting <http://astropy.readthedocs.org/en/latest/modeling/index.html>`_
to fit models to its spectra. For example, you can fit one model to the
continuum, another to an emission line of interest, and yet another to an
absorption line.

Currently, the following models are available:

========================= ======================================================
Pyfocal Model Name        Astropy Model Class
========================= ======================================================
BrokenPowerLaw            `~astropy.modeling.models.BrokenPowerLaw1D`
Const                     `~astropy.modeling.models.Const1D`
ExponentialCutoffPowerLaw `~astropy.modeling.models.ExponentialCutoffPowerLaw1D`
Gaussian                  `~astropy.modeling.models.Gaussian1D`
GaussianAbsorption        `~astropy.modeling.models.GaussianAbsorption1D`
Linear                    `~astropy.modeling.models.Linear1D`
LogParabola               `~astropy.modeling.models.LogParabola1D`
Lorentz                   `~astropy.modeling.models.Lorentz1D`
MexicanHat                `~astropy.modeling.models.MexicanHat1D`
Trapezoid                 `~astropy.modeling.models.Trapezoid1D`
PowerLaw                  `~astropy.modeling.models.PowerLaw1D`
Redshift                  `~astropy.modeling.models.Redshift`
Scale                     `~astropy.modeling.models.Scale`
Shift                     `~astropy.modeling.models.Shift`
Sine                      `~astropy.modeling.models.Sine1D`
Voigt                     `~astropy.modeling.models.Voigt1D`
========================= ======================================================

The models can be fitted with the following fitters:

=================== ============================================
Pyfocal Fitter Name Astropy Fitter Class
=================== ============================================
Levenberg-Marquardt `~astropy.modeling.fitting.LevMarLSQFitter`
Simplex             `~astropy.modeling.fitting.SimplexLSQFitter`
SLSQP               `~astropy.modeling.fitting.SLSQPLSQFitter`
=================== ============================================

To add a model:

#. Select the desired layer from "Layers" (left panel). For example, you can
   choose the layer containing your emission or absorption line.
   See :ref:`doc_viewer` on how to create a layer for ROI.
#. Select the desired model name from "Add Model" drop-down box and click "Add".
#. Scroll down (if needed) and click "Create".
#. A new model layer will be created under "Layers" (left panel).

To fine-tune model parameters:

#. Select the model layer under "Layers" (left panel) that contains the desired
   model.
#. Expand the model listing under "Current Models".
#. Double-click on the desired model parameter value in the listing.
#. When you see a blinking cursor, enter the new value and press "Enter".
#. Scroll down (if needed) and click "Update".

To fit a model:

#. Select the desired fitting routine from its drop-down menu.
#. Click "Perform Fit".

.. note::

    "Arithmetic Behavior" is used to define the relationship between different
    models for the same layer. It will be fully functional in a future release.
