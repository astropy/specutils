from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# STDLIB
import os
import logging

# LOCAL
from ..third_party.qtpy.QtCore import *
from ..third_party.qtpy.QtGui import *
from ..interfaces.managers import (data_manager, layer_manager,
                                   model_manager, plot_manager)
from ..interfaces.registries import loader_registry
from ..analysis.statistics import stats
from ..core.comms import Dispatch, DispatchHandle


class Controller(object):
    """
    GUI controller. This module only maintains the connections between
    gui elements.
    """
    def __init__(self, viewer):
        self.viewer = viewer

        # self._setup_events()
        self._setup_connections()
        self._setup_communications()

        # Initially, disable fitting controls so as to avoid forbidden states.
        # These buttons will be re-enabled on demand by the Viewer.
        self.viewer.main_window.comboBox_2.setEnabled(False)
        self.viewer.main_window.pushButton_3.setEnabled(False)

        DispatchHandle.setup(self)

    def _setup_communications(self):
        # Listen for subwindow selection events, update layer list on selection
        self.viewer.main_window.mdiArea.subWindowActivated.connect(
            lambda sw: Dispatch.on_select_window.emit(
                window=sw))

        # Listen for layer selection events, update model tree on selection
        self.viewer.wgt_layer_list.itemSelectionChanged.connect(
            lambda: Dispatch.on_select_layer.emit(
                layer_item=self.viewer.current_layer_item))

        # When a layer is selected, make that line more obvious than the others
        self.viewer.wgt_layer_list.itemSelectionChanged.connect(
            lambda: Dispatch.on_select_plot.emit(
                layer=self.viewer.current_layer))

        # When an interactable widget inside a layer item is clicked
        self.viewer.wgt_layer_list.itemChanged.connect(
            lambda li, col: Dispatch.on_clicked_layer.emit(
                layer_item=li))

        # Create a new layer based on any active ROIs
        self.viewer.main_window.toolButton_6.clicked.connect(
            lambda: self.add_roi_layer(self.viewer.current_layer,
                                       self.get_roi_mask(),
                                       self.viewer.current_sub_window))

    def _setup_connections(self):
        self.viewer.main_window.toolButton_3.clicked.connect(
            lambda: self.add_sub_window(
                data=self.viewer.current_data,
                window=None))

        self.viewer.main_window.actionOpen.triggered.connect(self.open_file)

        # When the arithmetic button is clicked, show math dialog
        self.viewer.main_window.arithmeticToolButton.clicked.connect(
            self._show_arithmetic_dialog)

        # Populate model dropdown
        self.viewer.main_window.comboBox.addItems(model_manager.all_models)

        # Populate fitting algorithm dropdown
        self.viewer.main_window.comboBox_2.addItems(model_manager.all_fitters)

        # When the add new model button is clicked, create a new model
        self.viewer.main_window.pushButton.clicked.connect(
            self.add_model)

        # Attach the fit button
        self.viewer.main_window.pushButton_3.clicked.connect(
            self.fit_model_layer)

        # Attach the create button
        self.viewer.main_window.pushButton_4.clicked.connect(
            self.add_model_layer)

        # Attach the update button
        self.viewer.main_window.pushButton_2.clicked.connect(
            self.update_model_layer)

    def _set_active_plot(self):
        current_sub_window = self.viewer.current_sub_window

        if current_sub_window is not None:
            current_sub_window.set_active_plot(
                self.viewer.current_layer)

    def _show_arithmetic_dialog(self):
        if self.viewer._layer_arithmetic_dialog.exec_():
            formula = self.viewer._layer_arithmetic_dialog\
                .ui_layer_arithmetic_dialog.formulaLineEdit.text()

            new_layer = layer_manager.add_from_formula(formula)

    def read_file(self, file_name):
        """
        Convenience method that directly reads a spectrum from a file.
        This exists mostly to facilitate development workflow. In time it
        could be augmented to support fancier features such as wildcards,
        file lists, mixed file types, and the like.
        Note that the filter string is hard coded here; its details might
        depend on the intrincacies of the registries, loaders, and data
        classes. In other words, this is brittle code.
        """
        file_name = str(file_name)
        file_ext = os.path.splitext(file_name)[-1]

        if file_ext in ('.txt', '.dat'):
            file_filter = 'ASCII (*.txt *.dat)'
        else:
            file_filter = 'Generic Fits (*.fits *.mits)'

        data = data_manager.load(file_name, file_filter)
        # self.viewer.add_data_item(data)

    def open_file(self, file_name=None):
        """
        Creates a `pyfocal.core.data.Data` object from the `Qt` open file
        dialog, and adds it to the data item list in the UI.
        """
        if not file_name or file_name is None:
            file_name, selected_filter = self.viewer.open_file_dialog(
                loader_registry.filters)

        if not file_name or file_name is None:
            return

        data = data_manager.load(str(file_name), str(selected_filter))

    def add_roi_layer(self, layer, mask=None, window=None):
        """
        Creates a layer object from the current ROIs of the active plot layer.

        Parameters
        ----------
        layer : pyfocal.core.data.Layer
            The current active layer of the active plot.
        window : QtGui.QMdiSubWindow
            The parent object within which the plot window resides.
        mask : ndarray
            Boolean mask.
        """
        roi_mask = mask if mask is not None else self.get_roi_mask(layer=layer)
        layer = layer_manager.new(layer._source,
                                  mask=roi_mask,
                                  window=window,
                                  name=layer._source.name + " Layer Slice")

        plot_container = plot_manager.new(layer, window)

    def add_sub_window(self, data=None, window=None, *args, **kwargs):
        """
        Creates a new plot widget to display in the MDI area. `data` and
        `sub_window` will be retrieved from the viewer if they are not defined.
        """
        data = data if data is not None else self.viewer.current_data

        if data is None:
            logging.warning("No data could be found with which to create new "
                            "sub window.")
            return

        if window is None:
            window = self.viewer.add_sub_window()

        layer = layer_manager.new(data, window=window)
        plot_container = plot_manager.new(layer, window)

    def add_model_layer(self):
        """
        Creates a new layer object using the currently defined model.
        """
        current_layer = self.viewer.current_layer
        model_inputs = self.viewer.get_model_inputs()

        # Remove models attached to parent layer
        model_manager.remove(current_layer)

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
            window=current_layer._window,
            name="New Model Layer",
            model=compound_model)

        plot_container = plot_manager.new(new_model_layer,
                                          current_layer._window)

        return new_model_layer

    def add_model(self):
        model_name = self.viewer.current_model
        layer = self.viewer.current_layer

        model = model_manager.new(model_name, layer)

        return model

    def fit_model_layer(self):
        current_layer = self.viewer.current_layer

        # Update the model parameters with those in the gui
        self.update_model_layer()

        # Create fitted layer
        fitted_layer = model_manager.fit_model(
            layer=current_layer,
            fitter_name=self.viewer.current_fitter)

    def get_roi_mask(self, layer=None):
        """
        Retrieves the array mask depending on the ROIs currently in the
        active plot window.

        Parameters
        ----------
        layer : Layer
            The layer containing the data from which the mask will be
            constructed.

        Returns
        -------
        roi_mask : ndarray
            A boolean array the size of currently selected layer masking
            outside the bounds of the ROIs.
        """
        current_layer = layer or self.viewer.current_layer
        current_sub_window = self.viewer.current_sub_window

        if current_sub_window is not None:
            roi_mask = current_sub_window.get_roi_mask(layer=current_layer)

            return roi_mask

    def update_model_layer(self, *args, **kwargs):
        """
        Updates the current layer with the results of the model.
        """
        current_layer = self.viewer.current_layer
        current_window = self.viewer.current_sub_window
        model_inputs = self.viewer.get_model_inputs()

        # Update model mask, only if rois exist
        mask = self.get_roi_mask(layer=current_layer._parent)
        mask = mask if len(current_window._rois) > 0 else None

        model_manager.update_model(layer=current_layer,
                                   model_inputs=model_inputs,
                                   formula=self.viewer.current_model_formula,
                                   mask=mask)

        # plot_manager.update_plots(current_window, current_layer)
        Dispatch.on_update_plot.emit(current_layer)

    @DispatchHandle.register_listener("on_clicked_layer")
    def _set_layer_visibility(self, layer_item, col=0):
        """
        Toggles the visibility of the plot in the sub window.

        Parameters
        ----------
        layer : Layer
            Layer object to toggle visibility.

        col : int
            QtTreeWidget data column.
        """
        layer = layer_item.data(0, Qt.UserRole)

        if layer is not None:
            current_window = layer._window

            current_window.set_visibility(
                layer, layer_item.checkState(col) == Qt.Checked, override=True)

    @DispatchHandle.register_listener("on_select_layer")
    def _update_layer_name(self, layer_item, col=0):
        layer = layer_item.data(0, Qt.UserRole)

        if hasattr(layer, 'name'):
            layer.name = layer_item.text(0)

        # Alert the statistics container to update the displayed layer name
        Dispatch.on_update_roi.emit(roi=None)

    @DispatchHandle.register_listener("on_select_model")
    def _update_model_name(self, model_item, col=0):
        model = model_item.data(0, Qt.UserRole)

        if hasattr(model, '_name'):
            model._name = model_item.text(0)

    @DispatchHandle.register_listener("on_select_layer", "on_update_roi")
    def update_statistics(self, layer_item=None, roi=None):
        if layer_item is not None:
            current_layer = layer_item.data(0, Qt.UserRole)
        else:
            current_layer = self.viewer.current_layer

        if current_layer is None:
            return

        mask = self.get_roi_mask(layer=current_layer)

        values = current_layer.data[mask]
        stat_dict = stats(values)

        Dispatch.on_update_stats.emit(stats=stat_dict, layer=current_layer)

    @DispatchHandle.register_listener("on_select_window")
    def update_layer_list(self, *args, **kwargs):
        """
        Clears and repopulates the layer list depending on the currently
        selected sub window.
        """
        current_window = self.viewer.current_sub_window
        layers = layer_manager.get_window_layers(current_window)
        self.viewer.clear_layer_widget()

        for layer in layers:
            container = plot_manager.get_plot_from_layer(layer=layer,
                                                         window=current_window)
            pixmap = QPixmap(10, 10)
            pixmap.fill(container._pen_stash['pen_on'].color())
            icon = QIcon(pixmap)
            self.viewer.add_layer_item(layer, icon=icon)

    @DispatchHandle.register_listener("on_select_layer", "on_update_model")
    def update_model_list(self, layer_item=None, model=None, layer=None):
        """
        Clears and repopulates the model list depending on the currently
        selected layer list item.
        """
        if layer_item is not None or layer is not None:
            current_layer = layer or self.viewer.current_layer
            self.viewer.clear_model_widget()

            if not hasattr(current_layer, 'model'):
                return

            if hasattr(current_layer.model, "submodel_names"):
                for i in range(len(current_layer.model.submodel_names)):
                    self.viewer.add_model_item(current_layer.model[i],
                                               current_layer)
            else:
                self.viewer.add_model_item(current_layer.model,
                                           current_layer)
        elif model is not None:
            current_layer = layer

            self.viewer.update_model_item(model)
