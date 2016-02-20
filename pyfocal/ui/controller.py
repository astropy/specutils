from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# STDLIB
import os
import logging

# LOCAL
from ..third_party.qtpy.QtCore import *
from ..interfaces.managers import (data_manager, layer_manager,
                                   model_manager, plot_manager)
from ..interfaces.registries import loader_registry
from ..interfaces import model_io
from ..analysis.statistics import stats
from ..analysis.modeling import apply_model

# To memorize last visited directory.
_model_directory = os.environ["HOME"]


class Controller(object):
    """GUI controller."""
    def __init__(self, viewer):
        # Controller-specific events
        self.viewer = viewer

        # self._setup_events()
        self._setup_connections()
        self._setup_model_fitting()

    # def _setup_events(self):
    #     # Setup layer event calling
    #     layer_manager.on_add += self.viewer.add_layer_item
    #     layer_manager.on_remove += self.viewer.remove_layer_item
    #
    #     # Setup model event calling
    #     model_layer_manager.on_add += self.viewer.add_layer_item
    #     model_layer_manager.on_remove += self.viewer.remove_layer_item
    #     model_layer_manager.on_add_model += self.viewer.add_model_item
    #     model_layer_manager.on_remove_model += self.viewer.remove_model_item

    def _setup_connections(self):
        self.viewer.main_window.actionOpen.triggered.connect(self.open_file)
        self.viewer.main_window.toolButton_3.clicked.connect(
            lambda: self.new_plot_sub_window())

        # Listen for subwindow selection events, update layer list on selection
        self.viewer.main_window.mdiArea.subWindowActivated.connect(
            self.update_layer_list)

        self.viewer.main_window.mdiArea.subWindowActivated.connect(
            self.update_model_list)

        # Listen for layer selection events, update model tree on selection
        self.viewer.wgt_layer_list.itemSelectionChanged.connect(
            self.update_model_list)

        # When a layer is selected, update the statistics for current ROIs
        self.viewer.wgt_layer_list.itemSelectionChanged.connect(
            self.update_statistics)

        # When a layer is selected, make that line more obvious than the others
        self.viewer.wgt_layer_list.itemSelectionChanged.connect(
            self._set_active_plot)

        # Create a new layer based on any active ROIs
        self.viewer.main_window.toolButton_6.clicked.connect(
                lambda: self.add_roi_layer(self.viewer.current_layer(),
                                           self.get_roi_mask(),
                                           self.viewer.current_sub_window()))

        # When clicking a layer, update the model list to show particular
        # buttons depending on if the layer is a model layer
        self.viewer.wgt_layer_list.itemSelectionChanged.connect(
            self.viewer._set_model_tool_options)

        # When clicking the checkbox next to a layer item, toggle the
        # visibility of that layer
        self.viewer.wgt_layer_list.itemClicked.connect(
            self._set_layer_visibility)

        # When the arithmetic button is clicked, show math dialog
        self.viewer.main_window.arithmeticToolButton.clicked.connect(
            self._show_arithmetic_dialog)

        # If a layer item is edited, make sure to save the name
        self.viewer.wgt_layer_list.itemChanged.connect(
            self._update_layer_name)

        # Attach the model save/read buttons
        self.viewer.main_window.pushButton_5.clicked.connect(
            self.save_model)
        # self.viewer.main_window.pushButton_6.clicked.connect(
        #     self.read_model)

    def _setup_sub_window_connections(self):
        # When the user changes the top axis
        pass

    def _setup_model_fitting(self):
        # Populate model dropdown
        self.viewer.main_window.comboBox.addItems(model_manager.all_models)

        # Populate fitting algorithm dropdown
        self.viewer.main_window.comboBox_2.addItems(model_manager.all_fitters)

        # Attach the add button
        self.viewer.main_window.pushButton.clicked.connect(
            self.create_new_model)

        # Attach the fit button
        self.viewer.main_window.pushButton_3.clicked.connect(
            self.fit_model)

        # Attach the create button
        self.viewer.main_window.pushButton_4.clicked.connect(
            self.new_model_layer)

        # Attach the update button
        self.viewer.main_window.pushButton_2.clicked.connect(
            self.update_model_layer)

        # Initially, disable fitting controls so as to avoid forbidden states.
        # These buttons will be re-enabled on demand by the Viewer.
        self.viewer.main_window.comboBox_2.setEnabled(False)
        self.viewer.main_window.pushButton_3.setEnabled(False)
        self.viewer.main_window.pushButton_5.setEnabled(False)
        self.viewer.main_window.pushButton_6.setEnabled(True)

        # If the a model item is edited, make sure to save the name
        self.viewer.wgt_model_list.itemChanged.connect(
            self._update_model_name
        )


    def _set_active_plot(self):
        current_sub_window = self.viewer.current_sub_window()

        if current_sub_window is not None:
            current_sub_window.set_active_plot(
                self.viewer.current_layer())

    def create_new_model(self):
        model = model_manager.new_model(self.viewer.current_layer(),
                                        self.viewer.current_model)

        self.viewer.add_model_item(model)

    def fit_model(self, *args):
        # when fitting, the selected layer is a ModelLayer, thus
        # the data to be fitted resides in the parent
        current_layer = self.viewer.current_layer()
        parent_layer = current_layer._parent
        mask = self.get_roi_mask(parent_layer)

        if parent_layer is None or mask is None:
            return

        # fit
        flux = parent_layer.data[mask]
        dispersion = parent_layer.dispersion[mask]

        model_dict = self.viewer.get_model_inputs()

        fitter_name = self.viewer.current_fitter
        formula = self.viewer.current_model_formula
        model = model_manager.get_compound_model(model_dict, formula=formula)

        # If the number of parameters is greater than the number of data
        # points, bail
        if len(model.parameters) > flux.size:
            logging.warning("Unable to perform fit; number of parameters is "
                            "greater than the number of data points.")
            return

        fitted_model = apply_model(model, dispersion, flux, fitter_name=fitter_name)

        # add fit results to current model layer.
        current_window = self.viewer.current_sub_window()

        current_layer.model = fitted_model

        # update GUI with fit results
        self.update_model_list()

        # plot fitted model
        plot_manager.update_plots(current_window, current_layer)

        # after model is fitted and plotted, the statistics window should
        # be refreshed to be consistent with the currently selected layer,
        # that is, the layer with the model just fitted.
        self.update_statistics()

    def save_model(self):
        model_dict = self.viewer.get_model_inputs()
        formula = self.viewer.current_model_formula
        model = model_manager.get_compound_model(model_dict, formula=formula)

        global _model_directory
        model_io.saveModelToFile(self.viewer.main_window.mdiArea, model, _model_directory)

    def update_statistics(self, *args):
        current_layer = self.viewer.current_layer()
        mask = self.get_roi_mask()

        if current_layer is None or mask is None:
            return

        values = current_layer.data[mask]

        stat_dict = stats(values)

        self.viewer.update_statistics(stat_dict)

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
        self.viewer.add_data_item(data)

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
        self.viewer.add_data_item(data)

    def add_roi_layer(self, layer, mask=None, sub_window=None):
        """
        Creates a layer object from the current ROIs of the active plot layer.

        Parameters
        ----------
        layer : pyfocal.core.data.Layer
            The current active layer of the active plot.
        sub_window : QtGui.QMdiSubWindow
            The parent object within which the plot window resides.
        mask : ndarray
            Boolean mask.
        """
        roi_mask = mask if mask is not None else self.get_roi_mask()
        layer = layer_manager.new_layer(layer._source,
                                        mask=roi_mask,
                                        window=sub_window,
                                        name=layer._source.name + " Layer Slice")

        self.add_plot(layer=layer)

    def new_plot_sub_window(self, data=None, window=None):
        """
        Creates a new plot widget to display in the MDI area. `data` and
        `sub_window` will be retrieved from the viewer if they are not defined.
        """
        data = data if data is not None else self.viewer.current_data()

        if data is None:
            return

        if window is None:
            window = self.viewer.add_sub_window()

            # Connect the statistics to the roi interactions
            window._plot_widget.on_roi_update += self.update_statistics

            # Connect the top axis change events
            # window._dynamic_axis.

        self.add_plot(data=data, window=window)

    def add_plot(self, data=None, layer=None, window=None):
        """
        Not defining the `layer` parameter will result in a new `layer`
        being created from the `data` object, and being attached to the
        `sub_window` object.

        Parameters
        ----------
        data : pyfocal.core.data.Data, optional
            The `Data` object from which a new layer (and plot) will be
            created.
        layer : pyfocal.core.data.Layer, option
            The `Layer` object from which a new visible plot will be created.
        window : pyfocal.ui.widgets.plot_sub_window.PlotSubWindow, optional
            The sub window to which a new layer object will be attached if
            there is no `layer` parameter defined.

        """
        if layer is None:
            if window is None:
                raise AttributeError("Sub window not defined.")
            elif data is None:
                raise AttributeError("Data is not defined")

            layer = layer_manager.new_layer(data, window=window)
        else:
            if window is None:
                window = layer._window

        plot_container = plot_manager.new_line_plot(layer, window)
        window.add_container(plot_container)
        self.update_layer_list()

    def new_model_layer(self):
        """
        Creates a new layer object using the currently defined model.
        """
        current_layer = self.viewer.current_layer()
        model_inputs = self.viewer.get_model_inputs()

        # Remove models attached to parent layer
        model_manager.remove(current_layer)

        if current_layer is None or not model_inputs:
            return

        compound_model = model_manager.get_compound_model(
            model_dict=model_inputs,
            formula=self.viewer.current_model_formula)

        # Create new layer using current ROI masks, if they exist
        mask = self.get_roi_mask()

        new_model_layer = layer_manager.new_model_layer(
            model=compound_model,
            data=current_layer._source,
            mask=mask,
            parent=current_layer,
            window=current_layer._window,
            name="New Model Layer")

        model_manager.add(compound_model, new_model_layer)

        self.viewer.add_layer_item(new_model_layer)

        self.add_plot(layer=new_model_layer)

        self.update_model_list()
        self.update_layer_list()

        return new_model_layer

    def update_model_layer(self):
        """
        Updates the current layer with the results of the model.
        """
        current_layer = self.viewer.current_layer()
        current_window = self.viewer.current_sub_window()
        model_inputs = self.viewer.get_model_inputs()
        model_manager.update_model(current_layer, model_inputs,
                                   self.viewer.current_model_formula)

        plot_manager.update_plots(current_window, current_layer)

    def update_layer_list(self):
        """
        Clears and repopulates the layer list depending on the currently
        selected sub window.
        """
        current_window = self.viewer.current_sub_window()

        layers = layer_manager.get_window_layers(current_window)

        self.viewer.clear_layer_widget()

        for layer in layers:
            self.viewer.add_layer_item(layer)

    def update_model_list(self):
        """
        Clears and repopulates the model list depending on the currently
        selected layer list item.
        """
        current_layer = self.viewer.current_layer()
        # models = model_manager.get_model_layers(current_layer, True)

        self.viewer.clear_model_widget()

        if current_layer is None or not hasattr(current_layer, "model"):
            return

        if hasattr(current_layer.model, "submodel_names"):
            for i in range(len(current_layer.model.submodel_names)):
                self.viewer.add_model_item(current_layer.model[i])
        else:
            self.viewer.add_model_item(current_layer.model)

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
        current_layer = layer or self.viewer.current_layer()
        current_sub_window = self.viewer.current_sub_window()

        if current_sub_window is not None:
            roi_mask = current_sub_window.get_roi_mask(layer=current_layer)

            return roi_mask

    def _set_layer_visibility(self, layer_item, col=0):
        """Toggles the visibility of the plot in the sub window.

        Parameters
        ----------
        layer : Layer
            Layer object to toggle visibility.

        col : int
            QtTreeWidget data column.
        """
        current_plot_window = self.viewer.current_sub_window()
        layer = layer_item.data(0, Qt.UserRole)

        if current_plot_window is not None and layer is not None:
            current_plot_window.set_visibility(
                layer, layer_item.checkState(col) == Qt.Checked, override=True)

    def _update_layer_name(self, layer_item, col=0):
        layer = layer_item.data(0, Qt.UserRole)

        if hasattr(layer, 'name'):
            layer.name = layer_item.text(0)

    def _update_model_name(self, model_item, col=0):
        model = model_item.data(0, Qt.UserRole)

        if hasattr(model, '_name'):
            model._name = model_item.text(0)

    def _show_arithmetic_dialog(self):
        current_window = self.viewer.current_sub_window()
        layers = layer_manager.get_window_layers(current_window)

        self.viewer._layer_arithmetic_dialog.clear_layers()
        self.viewer._layer_arithmetic_dialog.populate_layers(layers)

        if self.viewer._layer_arithmetic_dialog.exec_():
            ind1 = self.viewer._layer_arithmetic_dialog.first_layer
            ind2 = self.viewer._layer_arithmetic_dialog.second_layer
            op = self.viewer._layer_arithmetic_dialog.operator

            if op == 0:
                new_layer = layers[ind1] + layers[ind2]
            elif op == 1:
                new_layer = layers[ind1] - layers[ind2]
            elif op == 2:
                new_layer = layers[ind1] * layers[ind2]
            elif op == 3:
                new_layer = layers[ind1] / layers[ind2]

            layer_manager.add_layer(new_layer)
            self.add_plot(layer=new_layer)
