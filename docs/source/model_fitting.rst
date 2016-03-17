.. _doc_model_fitting:

Model Fitting
=============

SpecViz utilizes
`Astropy Models and Fitting <http://astropy.readthedocs.org/en/latest/modeling/index.html>`_
to fit models to its spectra. For example, you can fit one model to the
continuum, another to an emission line of interest, and yet another to an
absorption line.

Currently, the following models are available:

========================= ==========================================================
SpecViz Model Name        Astropy Model Class
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
Spline                    `~scipy.interpolate.UnivariateSpline`
Voigt                     `~astropy.modeling.functional_models.Voigt1D`
========================= ==========================================================

The models can be fitted with the following fitters:

=================== ============================================
SpecViz Fitter Name Astropy Fitter Class
=================== ============================================
Levenberg-Marquardt `~astropy.modeling.fitting.LevMarLSQFitter`
Simplex             `~astropy.modeling.fitting.SimplexLSQFitter`
SLSQP               `~astropy.modeling.fitting.SLSQPLSQFitter`
=================== ============================================

To add a model:

#. Select the desired layer from ``Layers`` (left panel). For example, you can
   choose the layer containing your emission or absorption line.
   See :ref:`doc_viewer` on how to create a layer for ROI.
#. Select the desired model name from the ``Add Model`` drop-down box and click
   ``Select`` to add it to ``Current Models``.
#. If desired, repeat the above step to add additional models.
#. Scroll down (if needed) and click ``Create Layer``.
#. A new model layer will be created under ``Layers`` (left panel) and it is
   attached to the selected data layer.

To fine-tune model parameters:

#. Select the model layer under ``Layers`` (left panel) that contains the desired
   model.
#. If desired, double-click on the model name to rename it. When you see a
   blinking cursor, enter its new name and press "Enter".
#. Expand the model listing under ``Current Models`` on the right of the viewer.
#. Double-click on the desired model parameter value in the listing.
   When you see a blinking cursor, enter the new value and press ``Enter``.
#. Scroll down (if needed) and click ``Update Layer``.

To fit a model:

#. Select the model layer under ``Layers`` that contains the model(s) you wish to
   fit to your data.
#. Select the desired fitter from ``Fitting Routine`` using its drop-down menu.
#. Click ``Perform Fit``. This may take up to a few seconds, depending on the
   complexity of the fit.
#. The associated model parameters will be adjusted accordingly.

The ``Arithmetic Behavior`` text box is used to define the relationship between
different models for the same layer. If nothing is defined, the default is to
add all the models together. To describe a non-default model relationship,
enter the model names and math operators, as shown in the examples below and
then press ``Create Layer`` or ``Update Layer`` to produce the compound model::

    Linear1 + Gaussian1

::

    Linear1 * Gaussian1

::

    Gaussian1 - Gaussian2

The entity that results from lumping together all the models, and combining them
either using the arithmetic behavior expression, or just adding them all together,
is called a "compound model".


Model names
^^^^^^^^^^^

When added to the ``Current Models`` list, a model will receive a default name
that is generated from the model type (as listed in the drop down model selector)
plus a running numerical suffix.

These names can be re-defined by clicking on the default name and typing a new
name. Note that re-defining names will require that any eventual expression in
the ``Arithmetic Behavior`` text box should be edited accordingly.

For now, we are limited to only alphanumeric characters (and no white spaces) when
re-naming models.


Spline model
^^^^^^^^^^^^

Note that the Spline model is of an intrinsically different nature than the
other models included in the drop down list of models. The Spline model, when
added to a pre-existing list of models, or when added by itself to an empty
list, will immediately be fitted to the data within the currently defined
Regions Of Interest. That is, being a linear model, there is no need to iterate
in search of a "best fit" spline. It is just computed once and for all, and kept
as part of the compound model that is built from the models in the list and the
arithmetic behavior expression.

This implies that, to change the regions of interest that define the spline,
one has no other way than removing the spline from the list of models. Then,
redefine the regions of interest, and add a new spline to the list. To change
a spline parameter, there is no need do discard the spline. Just do it in the
same way as with other models: just type in the new value for the parameter and
click on ``Update Layer``.

Subsequently, when the fitter iterates the compound model in search of a best
solution, the spline model will act as a constant. That is, it will be used to
compute the global result of the compound model, but its parameters won't be
accessed, and varied, by the fitter. Thus, the spline parameters are not fitted,
they are just a convenient mechanism that enables user access to the parameter's
values.

The documentation for the spline model cane be seen here:

http://docs.scipy.org/doc/scipy-0.16.0/reference/generated/scipy.interpolate.UnivariateSpline.html

Note that SpecViz provides access, at this point, to just two of the parameters
in the scipy implementation of th spline function. Pay special attention to the
``smooth`` parameter. SpecViz initializes it to a 'best guess' (``len(wavelength)``).
Too small of a value in here may cause the spline to enter an infinite loop.
Change the ``smooth`` value with care, trying to stay close to the default
value.


.. note::

    Model arithmetic is a work in progress.


Saving and Exporting Models to File
-----------------------------------

Selecting a model layer under "Layers" will enable the
:ref:`Save <doc_model_save>` (the "floppy disk" icon) and
:ref:`doc_model_export` (the "out the door" icon) buttons under
"Current Models" on the right of the viewer. Saving a model to a file will
enable SpecViz to read back that model into a new model layer. Exporting a model
to a file wil create a Python script in a ``.py`` file. This file can be
directly imported by Python in a command-line session.

Click on either button to get a file dialog window. Type in a file name.
If this file name does not end with the correct suffix, the suffix will
automatically be appended. Click "Save", or just the Return/Enter key.
The correct suffix for saved and exported files are ``.yaml`` and ``.py``,
respectively.


.. _doc_model_save:

Save and Load
^^^^^^^^^^^^^

Saving the model to a file works in the same way as :ref:`doc_model_export`.
The difference is that a saved model can be later read back into SpecViz via
the "Load" button (the "folder" icon), also under "Current Models".

For the "Load" button to be enabled, a data (spectrum) layer (not a model layer)
must be selected under "Layers" listing. The selected ``.yaml`` model file will
generate a model that will be attached to a new model layer associated under the
selected data layer.

The file is writen using the YAML format. Being a plain text file with a
self-explanatory structure, it can be edited at will by the user, e.g., to add
bounds, fixed flags, and ties to the model parameters. Note that these extra,
user-defined attributes, won't be accessible from SpecViz's user interface.
They will however, be accessible by the fitter when a fit is run on the
model. They will also be written ou correctly, either when saving or exporting
the model.

.. note::

    YAML format for saved models and usage of advanced features like bounds
    and fixed flags are work in progress.


.. _doc_model_export:

Export
^^^^^^

This will save the model in the currently selected model layer to a file
that can be directly imported by Python. The file is just a plain text
file with the model expressed recorded as a Python expression. The model
is associated to a variable named ``'model1'``.

The following example uses the ``'test3.py'`` file name, and a model comprised
of a constant and a gaussian:

.. code-block:: python

 >>> import test3
 >>> test3
 <module 'test3' from '/my/saved/models/test3.py'>
 >>> test3.model1
 <CompoundModel0(amplitude_0=0.297160787184, amplitude_1=2.25396100263, mean_1=15117.1710847, stddev_1=948.493577186)>
 >>> print(test3.model1)
 Model: CompoundModel0
 Inputs: ('x',)
 Outputs: ('y',)
 Model set size: 1
 Parameters:
      amplitude_0    amplitude_1      mean_1       stddev_1
     -------------- ------------- ------------- -------------
     0.297160787184 2.25396100263 15117.1710847 948.493577186

The file can be edited at will by the user, e.g., to add bounds, fixed flags,
and ties to the model parameters.

.. note::

    Security issues importing model this way into Python and usage of advanced
    features like bounds and fixed flags are work in progress.

