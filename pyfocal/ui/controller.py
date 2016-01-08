from qtpy.QtCore import *

from pyfocal.ui.widgets.plots.plot_window import PlotWindow
from ..interfaces.managers import data_manager, layer_manager, model_manager
from ..interfaces.registries import loader_registry


class Controller(object):
    def __init__(self, viewer):
        self.viewer = viewer

        # May change in the future: ask the controller to maintain a mapping of
        # sub windows and plots
        self.active_plots = {}

        self._setup_connections()
        self._setup_model_fitting()

    def _setup_connections(self):
        self.viewer.main_window.actionOpen.triggered.connect(self.open_file)
        self.viewer.main_window.toolButton_3.clicked.connect(
            self.create_plot_window)

        # Listen for subwindow selection events, update layer list on selection
        self.viewer.main_window.mdiArea.subWindowActivated.connect(
            self.update_layer_list)

        self.viewer.main_window.mdiArea.subWindowActivated.connect(
            self.update_model_list)

        # Listen for layer selection events, update model tree on selection
        self.viewer.wgt_layer_list.itemSelectionChanged.connect(
            self.update_model_list)

        self.viewer.main_window.toolButton_6.clicked.connect(
                lambda: self.add_roi_layer(self.viewer.current_layer(),
                                           self.get_roi_mask(),
                                           self.viewer.current_sub_window()))

    def _setup_model_fitting(self):
        # Populate model dropdown
        self.viewer.main_window.comboBox.addItems(model_manager.all_models)

        # Attach the add button
        self.viewer.main_window.pushButton.clicked.connect(
            self.create_new_model)

    def open_file(self):
        file_name, selected_filter = self.viewer.open_file_dialog(
            loader_registry.filters)

        data = data_manager.load(str(file_name), str(selected_filter))
        self.viewer.add_data_item(data)

    def create_sub_window(self):
        """
        Creates a new sub window with a graph object taken from the
        currently selected data list.
        """
        # Create sub window
        new_sub_window, wgt_sub_window = self.viewer.add_sub_window()

        return new_sub_window, wgt_sub_window

    def create_new_layer(self, data, mask=None, sub_window=None):
        """
        Creates a new layer from the selected data along with any current
        ROIs existing on the plot.
        """
        # Create the main layer for this sub window
        layer = layer_manager.new(data, mask=mask, sub_window=sub_window)
        self.viewer.add_layer_item(layer)
        self.update_layer_list()

        return layer

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
        layer = self.create_new_layer(layer._source, mask=roi_mask,
                                      sub_window=sub_window)
        self.add_plot(layer)

    def add_plot(self, layer):
        current_sub_window = self.viewer.current_sub_window()
        current_plot_window = self.active_plots[current_sub_window]

        current_plot_window.add_data(layer, visible=True)

    def create_plot_window(self):
        """
        Creates a new plot widget to display in the MDI area.
        """
        # Generate new sub window
        new_sub_window, wgt_sub_window = self.create_sub_window()

        # Grab data from list widget
        current_data = self.viewer.current_data()

        # Generate new data layer
        layer = self.create_new_layer(current_data,
                                      sub_window=new_sub_window)

        wgt_profile_plot = PlotWindow(layer, parent=wgt_sub_window)
        wgt_sub_window.gridLayout.addWidget(wgt_profile_plot)
        new_sub_window.show()

        # Archive the association of plot window and sub window
        self.active_plots[new_sub_window] = wgt_profile_plot

    def create_new_model(self):
        """
        Creates a new model for the selected layer.
        """
        # Create a new model for the layer
        current_layer = self.viewer.current_layer()
        current_model_name = self.viewer.main_window.comboBox.currentText()
        model = model_manager.add(current_layer, current_model_name)
        self.viewer.add_model_item(model, current_model_name)

    def update_layer_list(self):
        current_window = self.viewer.main_window.mdiArea.activeSubWindow()

        layers = layer_manager.get_sub_window_layers(current_window)

        self.viewer.clear_layer_widget()

        for layer in layers:
            self.viewer.add_layer_item(layer)

        if len(layers) > 0:
            self.viewer.wgt_layer_list.setCurrentRow(0)

    def update_model_list(self):
        current_layer = self.viewer.current_layer()
        models = model_manager.get_layer_models(current_layer)

        self.viewer.clear_model_widget()

        for model in models:
            self.viewer.add_model_item(model, model.__class__.__name__)

    def get_roi_mask(self):
        current_layer = self.viewer.current_layer()
        current_sub_window = self.viewer.current_sub_window()
        current_plot = self.active_plots[current_sub_window]

        roi_mask = current_plot.get_roi_mask(current_layer)

        return roi_mask


