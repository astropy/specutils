.. _doc_model_fitting:

Model Fitting
=============

Pyfocal utilizes
`Astropy Models and Fitting <http://astropy.readthedocs.org/en/latest/modeling/index.html>`_
to fit models to its spectra. For example, you can fit one model to the
continuum, another to an emission line of interest, and yet another to an
absorption line.

Currently, the following models are available:

========================= ==========================================================
Pyfocal Model Name        Astropy Model Class
========================= ==========================================================
BrokenPowerLaw            `~astropy.modeling.powerlaws.BrokenPowerLaw1D`
Const                     `~astropy.modeling.functional_models.Const1D`
ExponentialCutoffPowerLaw `~astropy.modeling.powerlaws.ExponentialCutoffPowerLaw1D`
Gaussian                  `~astropy.modeling.functional_models.Gaussian1D`
GaussianAbsorption        `~astropy.modeling.functional_models.GaussianAbsorption1D`
Linear                    `~astropy.modeling.functional_models.Linear1D`
LogParabola               `~astropy.modeling.powerlaws.LogParabola1D`
Lorentz                   `~astropy.modeling.functional_models.Lorentz1D`
MexicanHat                `~astropy.modeling.functional_models.MexicanHat1D`
Trapezoid                 `~astropy.modeling.functional_models.Trapezoid1D`
PowerLaw                  `~astropy.modeling.powerlaws.PowerLaw1D`
Redshift                  `~astropy.modeling.functional_models.Redshift`
Scale                     `~astropy.modeling.functional_models.Scale`
Shift                     `~astropy.modeling.functional_models.Shift`
Sine                      `~astropy.modeling.functional_models.Sine1D`
Voigt                     `~astropy.modeling.functional_models.Voigt1D`
========================= ==========================================================

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
#. If needed, double-click on the model name to rename it. When you see a blinking cursor, enter
   its new name and press "Enter".
#. Expand the model listing under "Current Models".
#. Double-click on the desired model parameter value in the listing. When you see a blinking
   cursor, enter the new value and press "Enter".
#. Scroll down (if needed) and click "Update".

To fit a model:

#. Select the desired fitting routine from its drop-down menu.
#. Click "Perform Fit".

The "Arithmetic Behavior" text box is used to define the relationship between different models
for the same layer. If nothing is defined, the default is to add all the models together.
If you have multiple models of the same name, you must rename them such that they have unique
names before defining the relationship.
To describe a non-default model relationship, enter the model names and math operators, as
shown in the examples below and then press "Create" or "Update" to produce the compound model::

    Linear + Gaussian

::

    Linear * Gaussian

::

    Gaussian1 - Gaussian2
