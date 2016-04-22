from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# STDLIB
import logging

# Third party
from astropy.units import spectral_density, spectral

# LOCAL
from ..third_party.qtpy.QtCore import *
from ..third_party.qtpy.QtWidgets import *
from ..interfaces.managers import (data_manager, window_manager, layer_manager,
                                   model_manager, plot_manager)
from ..interfaces.registries import loader_registry
from ..analysis import statistics
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
        # Listen for layer selection events, update model tree on selection
        self.viewer.wgt_layer_list.itemSelectionChanged.connect(
            lambda: Dispatch.on_selected_layer.emit(
                layer_item=self.viewer.current_layer_item))

        # When a layer is selected, make that line more obvious than the others
        self.viewer.wgt_layer_list.itemSelectionChanged.connect(
            lambda: Dispatch.on_selected_plot.emit(
                layer=self.viewer.current_layer))

        # When an interactable widget inside a layer item is clicked
        self.viewer.wgt_layer_list.itemClicked.connect(
            lambda li, col: Dispatch.on_clicked_layer.emit(
                layer_item=li))

        # When an interactable widget inside a layer item is clicked
        self.viewer.wgt_layer_list.itemChanged.connect(
            lambda li, col: Dispatch.on_changed_layer.emit(
                layer_item=li))

        # When the model items in the model tree change
        self.viewer.wgt_model_list.itemChanged.connect(
            lambda mi, col: Dispatch.on_changed_model.emit(
                model_item=mi))

    def _setup_connections(self):
        # When the layer list delete button is pressed
        self.viewer.main_window.layerRemoveButton.clicked.connect(
            self.remove_layer)

        # When the model list delete button is pressed
        self.viewer.main_window.modelRemoveButton.clicked.connect(
            self.remove_model)

        # When the arithmetic button is clicked, show math dialog
        self.viewer.main_window.arithmeticToolButton.clicked.connect(
            self._show_arithmetic_dialog)

        # Create a new layer based on any active ROIs
        self.viewer.main_window.toolButton_6.clicked.connect(
            lambda: self.add_roi_layer())

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
            lambda: self._change_plot_color(container)
        )

        self.viewer.layer_context_menu.exec_(
            self.viewer.wgt_layer_list.viewport().mapToGlobal(point))

    def _change_plot_color(self, container):
        col = QColorDialog.getColor(container._pen_stash['pen_on'].color(),
                                    self.viewer.wgt_layer_list)

        if col.isValid():
            container.pen = col

            Dispatch.on_updated_plot.emit(container=container)
        else:
            logging.warning("Color is not valid.")

    def _set_active_plot(self):
        current_sub_window = self.viewer.current_sub_window

        if current_sub_window is not None:
            current_sub_window.set_active_plot(
                self.viewer.current_layer)

    def _show_arithmetic_dialog(self):
        if self.viewer.current_layer is None:
            return

        if self.viewer._layer_arithmetic_dialog.exec_():
            formula = self.viewer._layer_arithmetic_dialog\
                .line_edit_formula.text()

            current_window = self.viewer.current_sub_window
            current_layers = window_manager.get_layers(current_window)
            new_layer = layer_manager.add_from_formula(formula,
                                                       layers=current_layers)

            if new_layer is None:
                logging.warning("Formula not valid.")
                return

            # If units match, plot the resultant on the same sub window,
            # otherwise create a new sub window to plot the spectra
            data_units_equiv = new_layer.data.unit.is_equivalent(
                current_window._plot_units[1],
                equivalencies=spectral_density(new_layer.dispersion))
            disp_units_equiv = new_layer.dispersion.unit.is_equivalent(
                current_window._plot_units[0], equivalencies=spectral())

            if data_units_equiv and disp_units_equiv:
                self.add_sub_window(layer=new_layer, window=current_window)
            else:
                logging.info("{} not equivalent to {}.".format(
                    new_layer.data.unit, current_window._plot_units[1]))
                self.add_sub_window(layer=new_layer)

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

    @DispatchHandle.register_listener("on_file_read")
    def read_file(self, file_name, file_filter=None):
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

        if file_filter is None:
            if file_ext in ('.txt', '.dat'):
                file_filter = 'ASCII (*.txt *.dat)'
            else:
                file_filter = 'Generic Fits (*.fits *.mits)'

        try:
            data = data_manager.load(file_name, file_filter)
        except:
            logging.error("Incompatible loader for selected data.")

    @DispatchHandle.register_listener("on_file_open")
    def open_file(self, file_name=None):
        """
        Creates a `specviz.core.data.Data` object from the `Qt` open file
        dialog, and adds it to the data item list in the UI.
        """
        if file_name is None:
            file_name, selected_filter = self.viewer.open_file_dialog(
                loader_registry.filters)

            self.read_file(file_name, file_filter=selected_filter)

    def add_roi_layer(self, layer=None, mask=None, window=None):
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
        layer = layer if layer is not None else self.viewer.current_layer

        # User attempts to slice before opening a file
        if layer is None:
            return

        window = window if window is not None else self.viewer.current_sub_window
        roi_mask = mask if mask is not None else self.get_roi_mask(layer=layer)

        new_layer = layer_manager.new(layer._source,
                                      mask=roi_mask,
                                      name=layer._source.name + " Layer Slice")

        window_manager.add(new_layer, window)
        plot_container = plot_manager.new(new_layer, window)

    @DispatchHandle.register_listener("on_add_window")
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
        else:
            layer_manager.add(layer)

        window_manager.add(layer, window)
        plot_container = plot_manager.new(layer, window)
        self.update_statistics()

    @DispatchHandle.register_listener("on_add_to_window")
    def add_to_sub_window(self, data=None, window=None, layer=None):
        """
        Adds the selected data set to the currently active sub window.
        """
        data = data or self.viewer.current_data
        window = window or self.viewer.current_sub_window

        if layer is None:
            layer = layer_manager.new(data)
        else:
            layer_manager.add(layer)

        window_manager.add(layer, window)
        plot_container = plot_manager.new(layer, window)
        self.update_statistics()

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

    def remove_layer(self, layer=None):
        current_layer = layer or self.viewer.current_layer

        if current_layer is None:
            return

        layer_manager.remove(layer=current_layer)
        window_manager.remove(layer=current_layer)
        plot_manager.remove(layer=current_layer)
        model_manager.remove(layer=current_layer)

    @DispatchHandle.register_listener("on_remove_data")
    def remove_data(self, current_data):
        if current_data is None:
            return

        data_manager.remove(data=current_data)

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
            current_window = window_manager.get(layer)

            current_window.set_visibility(
                layer, layer_item.checkState(col) == Qt.Checked, override=True)

    @DispatchHandle.register_listener("on_selected_layer", "on_changed_layer")
    def _update_layer_name(self, layer_item, col=0):
        if layer_item is None:
            return

        layer = layer_item.data(0, Qt.UserRole)

        if hasattr(layer, 'name'):
            layer.name = layer_item.text(0)

        # Alert the statistics container to update the displayed layer name
        Dispatch.on_updated_roi.emit(roi=None)

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

    @DispatchHandle.register_listener("on_selected_layer", "on_updated_roi")
    def update_statistics(self, layer_item=None, roi=None, measure_rois=None):
        if layer_item is not None:
            current_layer = layer_item.data(0, Qt.UserRole)
        else:
            current_layer = self.viewer.current_layer

        if current_layer is None:
            return

        if measure_rois is not None:
            # Set the active tab to measured
            self.viewer.main_window.statsTabWidget.setCurrentIndex(1)

            cont1_mask = self.get_roi_mask(roi=measure_rois[0])
            cont1_data = current_layer.data[cont1_mask]
            cont1_stat_dict = statistics.stats(cont1_data)

            cont2_mask = self.get_roi_mask(roi=measure_rois[2])
            cont2_data = current_layer.data[cont2_mask]
            cont2_stat_dict = statistics.stats(cont2_data)

            line_mask = self.get_roi_mask(roi=measure_rois[1])

            line = layer_manager.copy(current_layer)

            ew, flux, avg_cont = statistics.eq_width(cont1_stat_dict,
                                                     cont2_stat_dict,
                                                     line,
                                                     mask=line_mask)
            cent = statistics.centroid(line - avg_cont, mask=line_mask)

            print(cent)

            stat_dict = {"eq_width": ew, "centroid": cent, "flux": flux,
                         "avg_cont": avg_cont}
        else:
            # Set the active tab to basic
            self.viewer.main_window.statsTabWidget.setCurrentIndex(0)

            mask = self.get_roi_mask(layer=current_layer)

            if mask is None:
                values = current_layer.data
            else:
                values = current_layer.data[mask[current_layer._mask]]

            stat_dict = statistics.stats(values)

        Dispatch.on_updated_stats.emit(stats=stat_dict, layer=current_layer)

    @DispatchHandle.register_listener("on_selected_window")
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

    @DispatchHandle.register_listener("on_selected_layer", "on_updated_model")
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
