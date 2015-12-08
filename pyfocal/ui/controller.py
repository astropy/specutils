from __future__ import absolute_import, division, print_function
from ..factories.registries import loader_factory


class Controller(object):
    def __init__(self, viewer):
        self._viewer = viewer

    def create_sub_window(self):

        new_sub_window = self._viewer.add_sub_window()

        self._setup_connections()

    def _setup_connections(self):
        self._viewer.main_window.actionOpen.triggered.connect(self.open_file)

    def open_file(self):
        file_name = self._viewer.open_file_dialog(loader_factory.filters)



