from ..ui.widgets.plugin import Plugin
from ..third_party.qtpy.QtWidgets import *
from ..third_party.qtpy.QtCore import *
from ..core.comms import Dispatch, DispatchHandle


class DataListPlugin(Plugin):
    name = "Data List"

    def setup_ui(self):
        self.layout_vertical.setContentsMargins(11, 11, 11, 11)

        self.list_widget_data_list = QListWidget(self)

        self.layout_vertical.addWidget(self.list_widget_data_list)

        self.layout_horizontal = QHBoxLayout()

        self.button_create_sub_window = QToolButton(self)
        self.button_add_to_sub_window = QToolButton(self)
        self.button_remove_data = QToolButton(self)

        self.layout_horizontal.addWidget(self.button_create_sub_window)
        self.layout_horizontal.addWidget(self.button_add_to_sub_window)
        self.layout_horizontal.addStretch()
        self.layout_horizontal.addWidget(self.button_remove_data)

        self.layout_vertical.addLayout(self.layout_horizontal)

    def setup_connections(self):
        # Connect the create new sub window button
        self.button_create_sub_window.clicked.connect(
            lambda: Dispatch.on_add_window.emit(data=self.current_data))

        # Connect the add to current plot window button
        self.button_add_to_sub_window.clicked.connect(
            lambda: Dispatch.on_add_to_window.emit(
                data=self.current_data))

        # When the layer list delete button is pressed
        self.button_remove_data.clicked.connect(
            lambda: Dispatch.on_remove_data.emit(
            self.current_data))

    @property
    def current_data(self):
        """
        Returns the currently selected data object from the data list widget.

        Returns
        -------
        data : specviz.core.data.Data
            The `Data` object of the currently selected row.
        """
        data_item = self.list_widget_data_list.currentItem()

        if data_item is not None:
            data = data_item.data(Qt.UserRole)
            return data

    @DispatchHandle.register_listener("on_added_data")
    def add_data_item(self, data):
        """
        Adds a `Data` object to the loaded data list widget.

        Parameters
        ----------
        data : specviz.core.data.Data
            The `Data` object to add to the list widget.
        """
        new_item = QListWidgetItem(data.name, self.list_widget_data_list)
        new_item.setFlags(new_item.flags() | Qt.ItemIsEditable)

        new_item.setData(Qt.UserRole, data)

        self.list_widget_data_list.setCurrentItem(new_item)

    @DispatchHandle.register_listener("on_removed_data")
    def remove_data_item(self, data):
        data_item = self.get_data_item(data)

        self.list_widget_data_list.takeItem(
            self.list_widget_data_list.row(data_item))

    def get_data_item(self, data):
        for i in range(self.list_widget_data_list.count()):
            data_item = self.list_widget_data_list.item(0)

            if data_item.data(Qt.UserRole) == data:
                return data_item
