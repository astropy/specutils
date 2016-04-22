from abc import ABCMeta, abstractmethod, abstractproperty

from ...third_party.qtpy.QtWidgets import *
from ...third_party.qtpy.QtCore import *
from ...core.comms import Dispatch, DispatchHandle


class PluginMeta(type):
    pass


class Plugin(QDockWidget):
    """
    Base object for plugin infrastructure.
    """
    __meta__ = ABCMeta

    _all_plugins = []


    def __init__(self, parent=None):
        super(Plugin, self).__init__(parent)

        self._all_plugins.append(self)
        DispatchHandle.setup(self)

        # GUI Setup
        self.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)

        self.contents = QWidget()
        self.layout_vertical = QVBoxLayout(self.contents)

        self.setWidget(self.contents)

        self.setWindowTitle(self.name)
        self.setup_ui()
        self.setup_connections()

    def _set_name(self, value):
        if isinstance(value, str):
            self.name = value
        else:
            raise TypeError("Inappropriate type for 'name' property.")

    def _get_name(self):
        return self.name

    name = abstractproperty(_set_name, _get_name)

    @abstractmethod
    def setup_ui(self):
        raise NotImplementedError()

    @abstractmethod
    def setup_connections(self):
        raise NotImplementedError()

    def tool_button(cls, text="My tool button", on_click=None):
        new_tool_button = QToolButton(cls)
        new_tool_button.setText(text)

        if on_click is not None:
            new_tool_button.clicked.connect(on_click)

        cls.layout_vertical.addWidget(new_tool_button)

        return new_tool_button

    @staticmethod
    def push_button(self, text="My tool button", on_click=None):
        pass

    def list_widget(self, on_select=None):
        pass



