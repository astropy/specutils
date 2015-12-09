from __future__ import absolute_import, division, print_function

from ..interfaces.registries import loader_registry
from ..interfaces.managers import data_manager
from .widgets.profile import Profile


class Controller(object):
    def __init__(self, viewer):
        self._viewer = viewer
        self._setup_connections()

    def _setup_connections(self):
        self._viewer.main_window.actionOpen.triggered.connect(self.open_file)
        self._viewer.main_window.actionNew_Window.triggered.connect(self.create_sub_window)

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
        new_sub_window = self._viewer.add_sub_window()
        current_data = self._viewer.current_data()
        # wgt_profile_plot = Profile(current_data)
        # new_sub_windowsetWidgettwgt_profile_plot()




