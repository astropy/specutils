from abc import ABCMeta, abstractmethod, abstractproperty

from ...third_party.qtpy.QtWidgets import *
from ...third_party.qtpy.QtCore import *
from ...third_party.qtpy.QtGui import *
from ...core.comms import Dispatch, DispatchHandle

from ...ui.widgets.utils import ICON_PATH


class PluginMeta(type):
    pass


class Plugin(QDockWidget):
    """
    Base object for plugin infrastructure.
    """
    __meta__ = ABCMeta

    location = 'hidden'
    priority = 1

    def __init__(self, parent=None):
        super(Plugin, self).__init__(parent)
        # Initialize this plugin's actions list
        self._actions = []

        # Keep a reference to the active sub window
        self._active_window = None
        self._current_layer = None

        DispatchHandle.setup(self)

        # GUI Setup
        self.setAllowedAreas(Qt.AllDockWidgetAreas)

        self.scroll_area = QScrollArea(self)
        self.scroll_area.setFrameShape(QFrame.NoFrame)
        self.scroll_area.setFrameShadow(QFrame.Plain)
        self.scroll_area.setLineWidth(0)
        self.scroll_area.setWidgetResizable(True)
        self.scroll_area.setGeometry(QRect(0, 0, 306, 553))

        # The main widget inside the scroll area
        self.contents = QWidget()
        self.layout_vertical = QVBoxLayout(self.contents)
        self.layout_vertical.setContentsMargins(11, 11, 11, 11)
        self.layout_vertical.setSpacing(6)

        self.scroll_area.setWidget(self.contents)

        self.setWidget(self.scroll_area)

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

    def add_tool_bar_actions(self, icon_path, name="", category=None,
                             description="", priority=0, enabled=True,
                             callback=None):
        action = QAction(self)
        icon = QIcon(icon_path)
        action.setIcon(icon)
        action.setIconText(name)
        action.setStatusTip(description)
        action.setEnabled(enabled)

        action.triggered.connect(callback if callback is not None else
                                 lambda: None)

        self._actions.append(dict(action=action,
                                  category=(category, 0) if not isinstance(
                                      category, tuple) else category,
                                  priority=priority))

        return action

    @property
    def active_window(self):
        return self._active_window

    @property
    def current_layer(self):
        return self._current_layer

    @DispatchHandle.register_listener("on_activated_window")
    def set_active_window(self, window):
        self._active_window = window

    @DispatchHandle.register_listener("on_selected_layer")
    def set_active_layer(self, layer_item):
        if layer_item is not None:
            self._current_layer = layer_item.data(0, Qt.UserRole)
        else:
            self._current_layer = None
