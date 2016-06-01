from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# STDLIB
import logging

# Third party

# LOCAL
from ..third_party.qtpy.QtCore import *
from ..third_party.qtpy.QtWidgets import *
from ..core.comms import Dispatch, DispatchHandle

# We pick up the desired format for model files here.
# In a future release we may want to use both formats,
# the YAML format to Save and Load models, and the .py
# format to export models that can be directly imported
# by scripts elsewhere.

from ..interfaces.model_io import yaml_model_io
from ..interfaces.model_io import py_model_io


# To memorize last visited directory.
_model_directory = os.environ["HOME"]


class Controller(object):
    """
    GUI controller. This module only maintains the connections between
    gui elements.
    """
    def __init__(self, viewer):
        self.viewer = viewer

        # self._setup_events()
        # self._setup_connections()
        # self._setup_communications()
        # self._setup_context_menus()

        DispatchHandle.setup(self)

    def _setup_communications(self):

        # When the model items in the model tree change
        self.viewer.wgt_model_list.itemChanged.connect(
            lambda mi, col: Dispatch.on_changed_model.emit(
                model_item=mi))

    def _setup_connections(self):

        # When the model list delete button is pressed
        self.viewer.main_window.modelRemoveButton.clicked.connect(
            self.remove_model)

        # Populate model dropdown
        self.viewer.main_window.modelsComboBox.addItems(
            model_manager.all_models)

        # Populate fitting algorithm dropdown
        self.viewer.main_window.fittingRoutinesComboBox.addItems(
            model_manager.all_fitters)

        # When the add new model button is clicked, create a new model
        self.viewer.main_window.addModelButton.clicked.connect(
            self.add_model)

        # Attach the fit button
        self.viewer.main_window.fitModelLayerButton.clicked.connect(
            self.fit_model_layer)

        # Attach the create button
        self.viewer.main_window.createModelLayerButton.clicked.connect(
            self.add_model_layer)

        # Attach the update button
        self.viewer.main_window.updateModelLayerButton.clicked.connect(
            self.update_model_layer)

        # Attach the model save/read buttons
        self.viewer.main_window.saveModelButton.clicked.connect(
            self.save_model)

        self.viewer.main_window.loadModelButton.clicked.connect(
            self.load_model)

        self.viewer.main_window.exportModelButton.clicked.connect(
            self.export_model)

    def _setup_context_menus(self):
        self.viewer.wgt_layer_list.customContextMenuRequested.connect(
            self._layer_context_menu)

    def _layer_context_menu(self, point):
        layer_item = self.viewer.wgt_layer_list.itemAt(point)

        if layer_item is None:
            return

        container = plot_manager.get_plot_from_layer(
            layer_item.data(0, Qt.UserRole),
            self.viewer.current_sub_window)

        self.viewer.layer_context_menu.act_change_color.triggered.disconnect()
        self.viewer.layer_context_menu.act_change_color.triggered.connect(

        )

        self.viewer.layer_context_menu.exec_(
            self.viewer.wgt_layer_list.viewport().mapToGlobal(point))

    def _set_active_plot(self):
        current_sub_window = self.viewer.current_sub_window

        if current_sub_window is not None:
            current_sub_window.set_active_plot(
                self.viewer.current_layer)

    def _prepare_model_for_save(self):
        model_dict = self.viewer.get_model_inputs()
        formula = self.viewer.current_model_formula

        if len(model_dict) == 0:
            return None, None

        return model_manager.get_compound_model(model_dict, formula=formula), formula

    def save_model(self):
        model, formula = self._prepare_model_for_save()

        if model:
            global _model_directory
            yaml_model_io.saveModelToFile(self.viewer.main_window.mdiArea, model,
                                          _model_directory, expression=formula)

    def export_model(self):
        model, formula = self._prepare_model_for_save()

        if model:
            global _model_directory
            py_model_io.saveModelToFile(self.viewer.main_window.mdiArea, model,
                                          _model_directory, expression=formula)

    def load_model(self):
        global _model_directory
        fname = QFileDialog.getOpenFileNames(self.viewer.main_window.mdiArea,
                                             'Read model file',
                                             _model_directory,
                                              yaml_model_io.MODEL_FILE_FILTER)

        # File dialog returns a tuple with a list of file names.
        # We get the first name from the first tuple element.
        if len(fname[0]) < 1:
            return
        fname = fname[0][0]

        compound_model, formula, _model_directory = yaml_model_io.buildModelFromFile(fname)

        # Put new model in its own sub-layer under current layer.
        current_layer = self.viewer.current_layer

        if current_layer is None:
            return

        # Create new model layer using current ROI masks, if they exist
        mask = self.get_roi_mask(layer=current_layer)

        new_model_layer = layer_manager.new(
            data=current_layer._source,
            mask=mask,
            parent=current_layer,
            name="New Model Layer",
            model=compound_model)

        current_window = self.viewer.current_sub_window
        window_manager.add(new_model_layer, current_window)

        # Add the models from the just read compound model to the new model layer.
        # Note that a single model component requires a slight different handling
        # technique, since it is not iterable as a compound model is.
        if not hasattr(compound_model, '_format_expression'):
            model_manager.add(model=compound_model, layer=new_model_layer)
        else:
            for model in compound_model._submodels:
                model_manager.add(model=model, layer=new_model_layer)

        plot_container = plot_manager.new(new_model_layer,
                                          current_window)

        # put formula in text edit widget
        self.viewer.main_window.lineEdit.setText(formula)

    def add_model_layer(self):
        """
        Creates a new layer object using the currently defined model.
        """
        current_layer = self.viewer.current_layer
        current_window = self.viewer.current_sub_window
        model_inputs = self.viewer.get_model_inputs()

        if current_layer is None or not model_inputs:
            return

        compound_model = model_manager.get_compound_model(
            model_dict=model_inputs,
            formula=self.viewer.current_model_formula)

        # Create new layer using current ROI masks, if they exist
        mask = self.get_roi_mask(layer=current_layer)

        new_model_layer = layer_manager.new(
            data=current_layer._source,
            mask=mask,
            parent=current_layer,
            name="New Model Layer",
            model=compound_model)

        window_manager.add(new_model_layer, current_window)

        # Transfer models attached to parent layer
        model_manager.transfer_models(current_layer, new_model_layer)

        current_window = self.viewer.current_sub_window
        plot_container = plot_manager.new(new_model_layer,
                                          current_window)

        return new_model_layer

    def add_model(self):
        model_name = self.viewer.current_model
        layer = self.viewer.current_layer

        if layer is None:
            return

        if hasattr(layer, '_model'):
            mask = self.get_roi_mask(layer._parent)
        else:
            mask = self.get_roi_mask(layer)

        model = model_manager.new(model_name, layer, mask)

        return model

    def fit_model_layer(self):
        current_layer = self.viewer.current_layer

        # Update the model parameters with those in the gui
        # self.update_model_layer()

        # Create fitted layer
        fitted_layer = model_manager.fit_model(
            layer=current_layer,
            fitter_name=self.viewer.current_fitter)

        plot_manager.update_plots(layer=fitted_layer)

    def update_model_layer(self, model_item=None):
        """
        Updates the values in the model.
        """
        current_layer = self.viewer.current_layer
        current_window = self.viewer.current_sub_window
        model_inputs = self.viewer.get_model_inputs()

        if current_layer is None or current_window is None or not model_inputs:
            return

        # Update model mask, only if rois exist
        mask = self.get_roi_mask(layer=current_layer._parent) \
               if len(current_window._rois) > 0 else None

        model_manager.update_model(layer=current_layer,
                                   model_inputs=model_inputs,
                                   formula=self.viewer.current_model_formula,
                                   mask=mask)

        plot_manager.update_plots(layer=current_layer)

    def remove_model(self):
        model_item = self.viewer.current_model_item

        if model_item is None:
            return

        model = model_item.data(0, Qt.UserRole)
        layer = self.viewer.current_layer

        model_manager.remove(layer=layer, model=model)

    @DispatchHandle.register_listener("on_selected_layer", "on_changed_layer")
    def _update_layer_name(self, layer_item, col=0):
        if layer_item is None:
            return

        layer = layer_item.data(0, Qt.UserRole)

        if hasattr(layer, 'name'):
            layer.name = layer_item.text(0)

        # Alert the statistics container to update the displayed layer name
        Dispatch.on_updated_rois.emit(rois=None)

    @DispatchHandle.register_listener("on_selected_model", "on_changed_model")
    def _update_model_name(self, model_item, col=0):
        if model_item is None:
            return

        model = model_item.data(0, Qt.UserRole)

        if hasattr(model, '_name'):
            model._name = model_item.text(0)

    @DispatchHandle.register_listener("on_changed_model")
    def _update_model_parameters(self, *args, **kwargs):
        current_layer = self.viewer.current_layer
        model_inputs = self.viewer.get_model_inputs()

        if len(model_inputs) > 0:
            model_manager.update_model(current_layer, model_inputs)

    # @DispatchHandle.register_listener("on_selected_window")
    def update_layer_list(self, window=None, layer=None):
        """
        Clears and repopulates the layer list depending on the currently
        selected sub window.
        """
        if window is None:
            window = self.viewer.current_sub_window
        elif isinstance(window, QMdiSubWindow):
            window = window.widget()

        self.viewer.wgt_layer_list.clear()

        layers = window_manager.get_layers(window)

        for layer in layers:
            container = plot_manager.get_plot_from_layer(layer=layer,
                                                         window=window)
            self.viewer.add_layer_item(layer, unique=True)
            self.viewer.update_layer_item(container)

    # @DispatchHandle.register_listener("on_selected_layer", "on_updated_model")
    def update_model_list(self, layer_item=None, model=None, layer=None):
        """
        Clears and repopulates the model list depending on the currently
        selected layer list item.
        """
        if layer_item is not None or layer is not None:
            current_layer = layer or self.viewer.current_layer
            self.viewer.wgt_model_list.clear()

            models = model_manager.get_models(current_layer)

            for model in models:
                self.viewer.add_model_item(model=model,
                                           layer=current_layer, unique=True)
        elif model is not None:
            self.viewer.update_model_item(model)
