from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import importlib, inspect

from ..third_party.qtpy.QtCore import *
from ..third_party.qtpy.QtWidgets import *
from ..third_party.qtpy.QtGui import *

from ..core.comms import Dispatch
from .widgets.windows import MainWindow
from .widgets.plugin import Plugin


class Viewer(object):
    """
    The `Viewer` is the main construction area for all GUI widgets. This
    object does **not** control the interactions between the widgets,
    but only their creation and placement.
    """
    def __init__(self):
        # Instantiate main window object
        self._all_tool_bars = {}

        self.main_window = MainWindow()
        self.menu_docks = QMenu("Plugins")
        self.main_window.menu_bar.addMenu(self.menu_docks)

        self.menu_docks.addSeparator()

        # self.main_window.setDockNestingEnabled(True)

        # Load system and user plugins
        self.load_plugins()

        # Setup up top-level connections
        self._setup_connections()

    def load_plugins(self):
        instance_plugins = []

        for mod in os.listdir(os.path.abspath(os.path.join(
                __file__, '..', '..', 'plugins'))):

            mod = mod.split('.')[0]
            mod = importlib.import_module("specviz.plugins." + mod)
            cls_members = inspect.getmembers(
                mod, lambda member: inspect.isclass(member)
                                    and Plugin in member.__bases__)

            if len(cls_members) == 0:
                continue

            for _, cls_plugin in cls_members:
                instance_plugins.append(cls_plugin(self.main_window))

        for instance_plugin in sorted(instance_plugins,
                                      key=lambda x: x.priority):
            if instance_plugin.location != 'hidden':
                if instance_plugin.location == 'right':
                    location = Qt.RightDockWidgetArea
                elif instance_plugin.location == 'top':
                    location = Qt.TopDockWidgetArea
                else:
                    location = Qt.LeftDockWidgetArea

                self.main_window.addDockWidget(location, instance_plugin)
                instance_plugin.show()

                # Add this dock's visibility action to the menu bar
                self.menu_docks.addAction(
                    instance_plugin.toggleViewAction())

            for action in instance_plugin._actions:
                tool_bar = self._get_tool_bar(action['category'])
                tool_bar.addAction(action['action'])

    def _get_tool_bar(self, name):
        if name is None:
            name = "User Plugins"

        if name not in self._all_tool_bars:
            tool_bar = self.main_window.addToolBar(name)
            tool_bar.setToolButtonStyle(Qt.ToolButtonTextUnderIcon)
            tool_bar.show()

            self._all_tool_bars[name] = tool_bar

        return self._all_tool_bars[name]

    def _setup_connections(self):
        # Listen for subwindow selection events, update layer list on selection
        self.main_window.mdi_area.subWindowActivated.connect(
            lambda wi: Dispatch.on_selected_window.emit(
            window=wi.widget() if wi is not None else None))



