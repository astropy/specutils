from __future__ import absolute_import, division, print_function

from ..interfaces.registries import loader_registry
from ..interfaces.managers import data_manager, layer_manager
from .widgets.baseplot import BasePlot

from PyQt4.QtCore import *


class Controller(object):

    def __init__(self, viewer):
        self._viewer = viewer

        self._setup_connections()

    def _setup_connections(self):
        self._viewer.main_window.actionOpen.triggered.connect(self.open_file)
        self._viewer.main_window.toolButton_3.clicked.connect(
            self.create_sub_window)

    def open_file(self):
        file_name, selected_filter = self._viewer.open_file_dialog(
            loader_registry.filters)

        data = data_manager.load(str(file_name), str(selected_filter))
        self._viewer.add_data_item(data)

    def create_sub_window(self):
        """
        Creates a new sub window with a graph object taken from the
        currently selected data list.
        """
        # Create sub window
        new_sub_window, wgt_sub_window = self._viewer.add_sub_window()

        return new_sub_window, wgt_sub_window

    def create_new_layer(self, data, sub_window):
        """
        Creates a new layer from the selected data along with any current
        ROIs existing on the plot.
        """
        # Create the main layer for this sub window
        layer = layer_manager.new(data, sub_window)
        self._viewer.add_layer_item(layer)

        return layer

    def create_new_plot(self):
        # Generate new sub window
        new_sub_window, wgt_sub_window = self.create_sub_window()

        # Grab data from list widget
        current_data = self._viewer.current_data()

        # Generate new data layer
        self.create_new_layer(current_data, wgt_sub_window)

        wgt_profile_plot = BasePlot(current_data, parent=wgt_sub_window)
        wgt_sub_window.gridLayout.addWidget(wgt_profile_plot)
        new_sub_window.show()


