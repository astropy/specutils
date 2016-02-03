from qtpy.QtCore import *

from ..interfaces.managers import data_manager, layer_manager, \
    model_manager, plot_manager
from ..interfaces.registries import loader_registry
from ..analysis.statistics import stats


class Controller(object):
    def __init__(self, viewer):
        self.viewer = viewer

        self._setup_events()
        self._setup_connections()
        self._setup_model_fitting()

    def _setup_events(self):
        # Setup layer event calling
        layer_manager.on_add += self.viewer.add_layer_item
        layer_manager.on_remove += self.viewer.remove_layer_item

        # Setup model event calling
        model_manager.on_add += self.viewer.add_model_item
        model_manager.on_remove += self.viewer.remove_model_item

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

        # When a layer is selected, make that line more obvious than the others
        self.viewer.wgt_layer_list.itemSelectionChanged.connect(
            lambda: self.viewer.current_sub_window().set_active_plot(
                self.viewer.current_layer()
            ))

        # Create a new layer based on any active ROIs
        self.viewer.main_window.toolButton_6.clicked.connect(
                lambda: self.add_roi_layer(self.viewer.current_layer(),
                                           self.get_roi_mask(),
                                           self.viewer.current_sub_window()))

    def _setup_model_fitting(self):
        # Populate model dropdown
        self.viewer.main_window.comboBox.addItems(model_manager.all_models)

        # Attach the add button
        self.viewer.main_window.pushButton.clicked.connect(
            lambda: model_manager.add(self.viewer.current_layer(),
                                      self.viewer.current_model))

        # Attach the create button
        self.viewer.main_window.pushButton_4.clicked.connect(
            self.new_model_layer)

        # Attach the update button
        self.viewer.main_window.pushButton_2.clicked.connect(
            self.update_model_layer)

    def update_statistics(self):
        # Grab all available rois
        current_layer = self.viewer.current_layer()
        mask = self.get_roi_mask()

        if current_layer is None or mask is None:
            return

        stat_dict = stats(current_layer.data[mask])

        self.viewer.update_statistics(stat_dict)

    def open_file(self):
        """
        Creates a `pyfocal.core.data.Data` object from the `Qt` open file
        dialog, and adds it to the data item list in the UI.
        """
        file_name, selected_filter = self.viewer.open_file_dialog(
            loader_registry.filters)

        if file_name is None:
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
        layer = layer_manager.add(layer._source, mask=roi_mask,
                                  parent=sub_window)

        self.add_plot(layer=layer)

    def new_plot_sub_window(self, data=None, sub_window=None):
        """
        Creates a new plot widget to display in the MDI area. `data` and
        `sub_window` will be retrieved from the viewer if they are not defined.
        """
        data = data if data is not None else self.viewer.current_data()
        sub_window = sub_window if sub_window is not None else \
            self.viewer.add_sub_window()

        self.add_plot(data=data, sub_window=sub_window)

    def add_plot(self, data=None, layer=None, sub_window=None):
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
        sub_window : pyfocal.ui.widgets.plot_sub_window.PlotSubWindow, optional
            The sub window to which a new layer object will be attached if
            there is no `layer` parameter defined.

        """
        if layer is None:
            if sub_window is None:
                raise AttributeError("Sub window not defined.")
            elif data is None:
                raise AttributeError("Data is not defined")

            layer = layer_manager.add(data, parent=sub_window)
        else:
            if sub_window is None:
                sub_window = layer._parent

        plot_container = plot_manager.new_line_plot(layer, sub_window)
        sub_window.add_container(plot_container)

    def new_model_layer(self):
        """
        Creates a new layer object using the currently defined model.
        """
        current_layer = self.viewer.current_layer()

        if current_layer is None:
            return

        model_layer = layer_manager.new(current_layer._source,
                                        mask=self.get_roi_mask(),
                                        parent=current_layer._parent)

        model_inputs = self.viewer.get_model_inputs()
        compound_model = model_manager.get_compound_model(layer=current_layer,
                                                          model_dict=model_inputs,
                                                          formula=self.viewer.current_model_formula)

        new_layer = layer_manager.add_from_model(model_layer,
                                                 compound_model)
        self.add_plot(layer=new_layer)

        return new_layer

    def update_model_layer(self):
        """
        Updates the current layer with the results of the model.
        """
        current_layer = self.viewer.current_layer()
        model_inputs = model_manager.get_compound_model(current_layer)
        layer_manager.update_layer(current_layer, model_inputs)

    def update_layer_list(self):
        """
        Clears and repopulates the layer list depending on the currently
        selected sub window.
        """
        current_window = self.viewer.current_sub_window()

        layers = layer_manager.get_sub_window_layers(current_window)

        self.viewer.clear_layer_widget()

        for layer in layers:
            self.viewer.add_layer_item(layer)

        if len(layers) > 0:
            self.viewer.wgt_layer_list.setCurrentRow(0)

    def update_model_list(self):
        """
        Clears and repopulates the model list depending on the currently
        selected layer list item.
        """
        current_layer = self.viewer.current_layer()
        models = model_manager.get_layer_models(current_layer)

        self.viewer.clear_model_widget()

        for model in models:
            self.viewer.add_model_item(current_layer, model,
                                       model.__class__.__name__)

    def get_roi_mask(self):
        """
        Retrieves the array mask depending on the ROIs currently in the
        active plot window.

        Returns
        -------
        roi_mask : ndarray
            A boolean array the size of currently selected layer masking
            outside the bounds of the ROIs.
        """
        current_layer = self.viewer.current_layer()
        current_sub_window = self.viewer.current_sub_window()

        roi_mask = current_sub_window.get_roi_mask(layer=current_layer)

        return roi_mask

        self.main_window.actionChange_Color.triggered.connect(
                self._color_dialog)

    def _color_dialog(self):
        color = QColorDialog.getColor()

        if color.isValid():
            print(color.name())
