from ..ui.widgets.plugin import Plugin
from ..third_party.qtpy.QtWidgets import *
from ..third_party.qtpy.QtCore import *
from ..third_party.qtpy.QtGui import *
from ..core.comms import Dispatch, DispatchHandle
from ..ui.widgets.utils import ICON_PATH
from ..core.data import GenericSpectrum1D
from ..core.threads import FileLoadThread

import logging

import astropy.io.registry as io_registry


class DataListPlugin(Plugin):
    name = "Data List"
    location = "left"

    def __init__(self, *args, **kwargs):
        super(DataListPlugin, self).__init__(*args, **kwargs)

        self.file_load_thread = FileLoadThread()

        self.file_load_thread.status.connect(
            Dispatch.on_status_message.emit)

        self.file_load_thread.result.connect(
            self._data_loaded)

        # Add tool tray buttons
        self.button_open_data = self.add_tool_bar_actions(
            name="Open",
            description='Open data file',
            icon_path=os.path.join(ICON_PATH, "Open Folder-48.png"),
            category=('Loaders', 5),
            priority=1,
            callback=lambda: Dispatch.on_file_open.emit())

    def setup_ui(self):
        UiDataListPlugin(self)

    def setup_connections(self):
        # Enable/disable buttons depending on selection
        self.list_widget_data_list.itemSelectionChanged.connect(
            self.toggle_buttons)

        # Connect the create new sub window button
        self.button_create_sub_window.clicked.connect(
            lambda: Dispatch.on_add_window.emit(data=self.current_data))

        # Connect the add to current plot window button
        self.button_add_to_sub_window.clicked.connect(
            lambda: Dispatch.on_add_to_window.emit(
                data=self.current_data))

        # When the data list delete button is pressed
        self.button_remove_data.clicked.connect(
            lambda: self.remove_data_item())

    def _data_loaded(self, data):
        Dispatch.on_added_data.emit(data=data)

        if self.active_window is None:
            Dispatch.on_add_window.emit(data=data)

    @property
    def current_data(self):
        """
        Returns the currently selected data object from the data list widget.

        Returns
        -------
        data : specviz.core.data.GenericSpectrum1D
            The `Data` object of the currently selected row.
        """
        data_item = self.list_widget_data_list.currentItem()

        if data_item is not None:
            data = data_item.data(Qt.UserRole)
            return data

    @property
    def current_data_item(self):
        return self.list_widget_data_list.currentItem()

    @DispatchHandle.register_listener("on_file_open")
    def open_file(self, file_name=None):
        """
        Creates a `specviz.core.data.Data` object from the `Qt` open file
        dialog, and adds it to the data item list in the UI.
        """
        if file_name is None:
            file_name, selected_filter = self.open_file_dialog()

            self.read_file(file_name, file_filter=selected_filter)

    def open_file_dialog(self):
        """
        Given a list of filters, prompts the user to select an existing file
        and returns the file path and filter.

        Returns
        -------
        file_name : str
            Path to the selected file.
        selected_filter : str
            The chosen filter (this indicates which custom loader from the
            registry to use).
        """
        dialog = QFileDialog(self)
        dialog.setFileMode(QFileDialog.ExistingFile)
        dialog.setNameFilters([x + " (*)" for x in
                               io_registry.get_formats(GenericSpectrum1D)[
                                   'Format']])

        if dialog.exec_():
            file_names = dialog.selectedFiles()
            selected_filter = dialog.selectedNameFilter().replace(" (*)", "")

            return file_names[0], selected_filter

        return None, None

    @DispatchHandle.register_listener("on_file_read")
    def read_file(self, file_name, file_filter=None):
        self.file_load_thread(file_name=file_name, file_filter=file_filter)
        self.file_load_thread.start()

    @DispatchHandle.register_listener("on_added_data")
    def add_data_item(self, data):
        """
        Adds a `Data` object to the loaded data list widget.

        Parameters
        ----------
        data : specviz.core.data.GenericSpectrum1D
            The `Data` object to add to the list widget.
        """
        new_item = QListWidgetItem(data.name, self.list_widget_data_list)
        new_item.setFlags(new_item.flags() | Qt.ItemIsEditable)

        new_item.setData(Qt.UserRole, data)

        self.list_widget_data_list.setCurrentItem(new_item)

    @DispatchHandle.register_listener("on_remove_data")
    def remove_data_item(self, data=None):
        if data is None:
            data = self.current_data

        data_item = self.get_data_item(data)

        self.list_widget_data_list.takeItem(
            self.list_widget_data_list.row(data_item))

        Dispatch.on_removed_data.emit(data=self.current_data)

    @DispatchHandle.register_listener("on_remove_all_data")
    def remove_all_data(self):
        self.list_widget_data_list.clear()

    def get_data_item(self, data):
        for i in range(self.list_widget_data_list.count()):
            data_item = self.list_widget_data_list.item(i)

            if data_item.data(Qt.UserRole) == data:
                return data_item

    def toggle_buttons(self):
        if self.current_data_item is not None:
            self.label_unopened.hide()
            self.button_remove_data.setEnabled(True)
            self.button_create_sub_window.setEnabled(True)
            # self.button_add_to_sub_window.setEnabled(True)
        else:
            self.label_unopened.show()
            self.button_remove_data.setEnabled(False)
            self.button_create_sub_window.setEnabled(False)
            # self.button_add_to_sub_window.setEnabled(False)


class UiDataListPlugin:
    def __init__(self, plugin):
        plugin.layout_vertical.setContentsMargins(11, 11, 11, 11)

        # List widget for the data sets
        plugin.list_widget_data_list = QListWidget(plugin)

        # Label box to show when no data set has been loaded
        plugin.label_unopened = QLabel(plugin)
        plugin.label_unopened.setAlignment(Qt.AlignCenter | Qt.AlignHCenter)
        plugin.label_unopened.setText("Click the folder icon to open a data set")
        plugin.label_unopened.setWordWrap(True)
        plugin.label_unopened.setStyleSheet("""
        QLabel {
            color: #8a6d3b;
            background-color: #fcf8e3;
            padding: 10px;
            border: 1px solid #faebcc;
            border-radius: 4px;
        }""")

        plugin.layout_vertical.addWidget(plugin.label_unopened)
        plugin.layout_vertical.addWidget(plugin.list_widget_data_list)

        plugin.layout_horizontal = QHBoxLayout()

        plugin.button_create_sub_window = QToolButton(plugin)
        plugin.button_create_sub_window.setIcon(QIcon(os.path.join(
            ICON_PATH, "Open in Browser-50.png")))
        plugin.button_create_sub_window.setIconSize(QSize(25, 25))
        plugin.button_create_sub_window.setMaximumSize(QSize(35, 35))
        plugin.button_create_sub_window.setEnabled(False)

        plugin.button_add_to_sub_window = QToolButton(plugin)
        plugin.button_add_to_sub_window.setIcon(QIcon(os.path.join(
            ICON_PATH, "Change Theme-50.png")))
        plugin.button_add_to_sub_window.setIconSize(QSize(25, 25))
        plugin.button_add_to_sub_window.setMaximumSize(QSize(35, 35))
        plugin.button_add_to_sub_window.setEnabled(False)

        plugin.button_remove_data = QToolButton(plugin)
        plugin.button_remove_data.setIcon(QIcon(os.path.join(
            ICON_PATH, "Delete-48.png")))
        plugin.button_remove_data.setEnabled(False)
        plugin.button_remove_data.setIconSize(QSize(25, 25))
        plugin.button_remove_data.setMaximumSize(QSize(35, 35))

        plugin.layout_horizontal.addWidget(plugin.button_create_sub_window)
        plugin.layout_horizontal.addWidget(plugin.button_add_to_sub_window)
        plugin.layout_horizontal.addStretch()
        plugin.layout_horizontal.addWidget(plugin.button_remove_data)

        plugin.layout_vertical.addLayout(plugin.layout_horizontal)


