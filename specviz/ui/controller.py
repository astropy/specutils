from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# STDLIB
import os
import logging

# LOCAL
from ..third_party.qtpy.QtCore import *
from ..third_party.qtpy.QtGui import *
from ..third_party.qtpy.QtWidgets import *
from ..interfaces.managers import (data_manager, window_manager, layer_manager,
                                   model_manager, plot_manager)
from ..interfaces.registries import loader_registry
from ..analysis import statistics
from ..core.comms import Dispatch, DispatchHandle
from ..interfaces import model_io

from astropy.units.core import UnitConversionError

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
        self._setup_connections()
        self._setup_communications()
        self._setup_context_menus()

        DispatchHandle.setup(self)

    def _setup_communications(self):
        # Listen for subwindow selection events, update layer list on selection
        self.viewer.main_window.mdiArea.subWindowActivated.connect(
            lambda sw: Dispatch.on_select_window.emit(
                window=sw))

        # Listen for layer selection events, update model tree on selection
        self.viewer.wgt_layer_list.currentItemChanged.connect(
            lambda ci, pi: Dispatch.on_select_layer.emit(
                layer_item=ci))

        # When a layer is selected, make that line more obvious than the others
        self.viewer.wgt_layer_list.itemSelectionChanged.connect(
            lambda: Dispatch.on_select_plot.emit(
                layer=self.viewer.current_layer))

        # When an interactable widget inside a layer item is clicked
        self.viewer.wgt_layer_list.itemClicked.connect(
            lambda li, col: Dispatch.on_clicked_layer.emit(
                layer_item=li))

        # When an interactable widget inside a layer item is clicked
        self.viewer.wgt_model_list.currentItemChanged.connect(
            lambda mi, col: Dispatch.on_changed_model.emit(
                model_item=mi))

        # Create a new layer based on any active ROIs
        self.viewer.main_window.toolButton_6.clicked.connect(
            lambda: self.add_roi_layer(self.viewer.current_layer,
                                       self.get_roi_mask(),
                                       self.viewer.current_sub_window))

    def _setup_connections(self):
        self.viewer.main_window.toolButton_3.clicked.connect(
            lambda: self.add_sub_window(
                data=self.viewer.current_data))

        self.viewer.main_window.actionOpen.triggered.connect(self.open_file)

        # When the layer list delete button is pressed
        self.viewer.main_window.layerRemoveButton.clicked.connect(
            lambda: layer_manager.remove(self.viewer.current_layer))

        # When the arithmetic button is clicked, show math dialog
        self.viewer.main_window.arithmeticToolButton.clicked.connect(
            self._show_arithmetic_dialog)

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
            self.update_model_list)

        # Attach the model save/read buttons
        self.viewer.main_window.saveModelButton.clicked.connect(
            self.save_model)
        self.viewer.main_window.loadModelButton.clicked.connect(
            self.load_model)

    def _setup_context_menus(self):
        self.viewer.wgt_layer_list.customContextMenuRequested.connect(
            self._layer_context_menu)

        self.viewer.wgt_model_list.customContextMenuRequested.connect(
            self._model_context_menu)

    def _layer_context_menu(self, point):
        layer_item = self.viewer.wgt_layer_list.itemAt(point)
        container = plot_manager.get_plot_from_layer(
            layer_item.data(0, Qt.UserRole),
            self.viewer.current_sub_window)

        self.viewer.layer_context_menu.act_change_color.triggered.disconnect()
        self.viewer.layer_context_menu.act_change_color.triggered.connect(
            lambda: self._change_plot_color(container)
        )

        self.viewer.layer_context_menu.exec_(
            self.viewer.wgt_layer_list.viewport().mapToGlobal(point))

    def _model_context_menu(self, point):
        model_item = self.viewer.wgt_model_list.itemAt(point)
        model = model_item.data(0, Qt.UserRole)
        layer = self.viewer.current_layer

        self.viewer.model_context_menu.act_remove.triggered.disconnect()
        self.viewer.model_context_menu.act_remove.triggered.connect(
            lambda: model_manager.remove(layer=layer, model=model)
        )

        self.viewer.model_context_menu.exec_(
            self.viewer.wgt_model_list.viewport().mapToGlobal(point))

    def _change_plot_color(self, container):
        col = QColorDialog.getColor(container._pen_stash['pen_on'].color(),
                                    self.viewer.wgt_layer_list)

        if col.isValid():
            container.pen = col

            plot_manager.update_plots(container=container)

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
            current_window = self.viewer.current_sub_window

            # If units match, plot the resultant on the same sub window,
            # otherwise create a new sub window to plot the spectra
            if new_layer.data.unit.is_equivalent(current_window._plot_units[1]):
                self.add_sub_window(layer=new_layer, window=current_window)
            else:
                print("{} not equivalent to {}.".format(
                    new_layer.data.unit, current_window._plot_units[1]))
                self.add_sub_window(layer=new_layer)

    def save_model(self):
        model_dict = self.viewer.get_model_inputs()
        formula = self.viewer.current_model_formula
        model = model_manager.get_compound_model(model_dict, formula=formula)

        global _model_directory
        model_io.saveModelToFile(self.viewer.main_window.mdiArea, model,
                                 _model_directory)

    def load_model(self):
        global _model_directory
        fname = QFileDialog.getOpenFileNames(self.viewer.main_window.mdiArea,
                                             'Read model file',
                                             _model_directory,
                                             "Pyhton files (*.py)")

        # File dialog returns a tuple with a list of file names.
        # We get the first name from the first tuple element.
        fname = fname[0][0]

        compound_model, _model_directory = model_io.buildModelFromFile(fname)

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
            window=current_layer._window,
            name="New Model Layer",
            model=compound_model)

        # Add the models from the just read compound model to the new model layer.
        # Note that a single model component requires a slight different handling
        # technique, since it is not iterable as a compound model is.
        if not hasattr(compound_model, '_format_expression'):
            model_manager.add(model=compound_model, layer=new_model_layer)
        else:
            for model in compound_model:
                model_manager.add(model=model, layer=new_model_layer)

        plot_container = plot_manager.new(new_model_layer,
                                          current_layer._window)

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
        Creates a `specviz.core.data.Data` object from the `Qt` open file
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
        layer : specviz.core.data.Layer
            The current active layer of the active plot.
        window : QtGui.QMdiSubWindow
            The parent object within which the plot window resides.
        mask : ndarray
            Boolean mask.
        """
        roi_mask = mask if mask is not None else self.get_roi_mask(layer=layer)
        layer = layer_manager.new(layer._source,
                                  mask=roi_mask,
                                  name=layer._source.name + " Layer Slice")

        window_manager.add(layer, window)
        plot_container = plot_manager.new(layer, window)

    def add_sub_window(self, data=None, window=None, layer=None):
        """
        Creates a new plot widget to display in the MDI area. `data` and
        `sub_window` will be retrieved from the viewer if they are not defined.
        """
        if data is None:
            data = self.viewer.current_data

        if window is None:
            window = self.viewer.add_sub_window()

        if layer is None:
            layer = layer_manager.new(data)

        window_manager.add(layer, window)

    def add_model_layer(self):
        """
        Creates a new layer object using the currently defined model.
        """
        current_layer = self.viewer.current_layer
        current_window = self.viewer.current_sub_window
        model_inputs = self.viewer.get_model_inputs()

        if current_layer is None or not model_inputs:
            return

        # Remove models attached to parent layer
        model_manager.remove(current_layer)

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

        # Add the models to the new model layer
        for model in model_inputs:
            model_manager.add(model=model, layer=new_model_layer)

        current_window = self.viewer.current_sub_window
        plot_container = plot_manager.new(new_model_layer,
                                          current_window)

        return new_model_layer

    def add_model(self):
        model_name = self.viewer.current_model
        layer = self.viewer.current_layer
        print(layer)
        model = model_manager.new(model_name, layer)

        return model

    def fit_model_layer(self):
        current_layer = self.viewer.current_layer

        # Update the model parameters with those in the gui
        # self.update_model_layer()

        # Create fitted layer
        fitted_layer = model_manager.fit_model(
            layer=current_layer,
            fitter_name=self.viewer.current_fitter)

    def get_roi_mask(self, layer=None, roi=None):
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
            roi_mask = current_sub_window.get_roi_mask(layer=current_layer,
                                                       roi=roi)

            return roi_mask

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
        mask = self.get_roi_mask(layer=current_layer._parent)
        mask = mask if len(current_window._rois) > 0 else None

        model_manager.update_model(layer=current_layer,
                                   model_inputs=model_inputs,
                                   formula=self.viewer.current_model_formula,
                                   mask=mask)

        plot_manager.update_plots(layer=self.viewer.current_layer)

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
            current_window = self.viewer.current_sub_window

            current_window.set_visibility(
                layer, layer_item.checkState(col) == Qt.Checked, override=True)

    @DispatchHandle.register_listener("on_select_layer", "on_clicked_layer")
    def _update_layer_name(self, layer_item, col=0):
        if layer_item is None:
            return

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

    @DispatchHandle.register_listener("on_changed_model")
    def _update_model_parameters(self, *args, **kwargs):
        current_layer = self.viewer.current_layer
        model_inputs = self.viewer.get_model_inputs()

        if len(model_inputs) > 0:
            model_manager.update_model(current_layer, model_inputs)

    @DispatchHandle.register_listener("on_added_layer_to_window")
    def plot_layer(self, layer=None, window=None):
        plot_container = plot_manager.new(layer=layer, window=window)

    @DispatchHandle.register_listener("on_remove_layer_from_window")
    def remove_plot_layer(self, layer=None, window=None):
        plot_manager.remove(layer=layer, window=window)

    @DispatchHandle.register_listener("on_select_layer", "on_update_roi")
    def update_statistics(self, layer_item=None, roi=None, measured_rois=None):
        if layer_item is not None:
            current_layer = layer_item.data(0, Qt.UserRole)
        else:
            current_layer = self.viewer.current_layer

        if current_layer is None:
            return

        if measured_rois is not None:
            cont1_mask = self.get_roi_mask(roi=measured_rois[0])
            cont1_data = current_layer.data[cont1_mask]
            cont1_stat_dict = statistics.stats(cont1_data)

            cont2_mask = self.get_roi_mask(roi=measured_rois[2])
            cont2_data = current_layer.data[cont2_mask]
            cont2_stat_dict = statistics.stats(cont2_data)

            line_mask = self.get_roi_mask(roi=measured_rois[1])

            line = layer_manager.copy(current_layer)
            line._mask = line_mask

            ew = statistics.eq_width(cont1_stat_dict, cont2_stat_dict, line)[1]

            stat_dict = {"eq_width": ew}
        else:
            mask = self.get_roi_mask(layer=current_layer)

            values = current_layer.data[mask]
            stat_dict = statistics.stats(values)

        Dispatch.on_update_stats.emit(stats=stat_dict, layer=current_layer)

    @DispatchHandle.register_listener("on_select_window",
                                      "on_added_layer_to_window")
    def update_layer_list(self, window=None, layer=None):
        """
        Clears and repopulates the layer list depending on the currently
        selected sub window.
        """
        if window is None:
            window = self.viewer.current_sub_window
        elif isinstance(window, QMdiSubWindow):
            window = window.widget()

        layers = window_manager.get_layers(window)
        self.viewer.clear_layer_widget()

        for layer in layers:
            container = plot_manager.get_plot_from_layer(layer=layer,
                                                         window=window)
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

            models = model_manager.get_models(current_layer)

            for model in models:
                self.viewer.add_model_item(model=model,
                                           layer=current_layer)
        elif model is not None:
            self.viewer.update_model_item(model)
