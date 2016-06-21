from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import importlib, inspect, sys, os

from ..third_party.qtpy.QtCore import *
from ..third_party.qtpy.QtWidgets import *
from ..third_party.qtpy.QtGui import *

from .widgets.sub_windows import PlotSubWindow
from ..core.comms import Dispatch, DispatchHandle
from ..core.annotation import LineIDMarker

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
        self._all_categories = {}
        self._all_tool_bars = {}

        self.main_window = MainWindow()
        self.menu_docks = QMenu("Plugins")
        self.main_window.menu_bar.addMenu(self.menu_docks)

        self.menu_docks.addSeparator()

        self.main_window.setDockNestingEnabled(True)

        # Load system and user plugins
        self.load_plugins()

        # Setup toolbar
        self.setup_toolbar()

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

                # Add this dock's visibility action to the menu bar
                self.menu_docks.addAction(
                    instance_plugin.toggleViewAction())

    def setup_toolbar(self):
        for action_dict in Plugin._tool_bar_actions:
            tool_bar = self.get_tool_bar(action_dict['category'])
            tool_bar.addAction(action_dict['widget'])

        # for button in Plugin._tool_buttons:
        #     tool_bar = self.get_tool_bar(button['category'])
        #     cat = self.get_category(button['category'])
        #     cat['grid'].addWidget(button['widget'])
        #
        # for cat_name, cat_dict in self._all_categories.items():
        #     tool_bar = self._all_tool_bars[cat_name]
        #     tool_bar.addWidget(cat_dict['widget'])

    def get_tool_bar(self, name):
        if name is None:
            name = "User Plugins"

        if name not in self._all_tool_bars:
            tool_bar = self.main_window.addToolBar(name)
            tool_bar.setToolButtonStyle(Qt.ToolButtonTextUnderIcon)

            label = QLabel()
            label.setText(name.lower())
            label.setStyleSheet("""
                QLabel {
                    color: #999999;
                    font-variant: small-caps;
                    font-weight: bold;
                    font-size: 0.7em;
                    border-top: 1px solid #999999;
                }""")

            self._all_tool_bars[name] = tool_bar

        return self._all_tool_bars[name]

    def _setup_connections(self):
        # Listen for subwindow selection events, update layer list on selection
        self.main_window.mdi_area.subWindowActivated.connect(
            lambda wi: Dispatch.on_selected_window.emit(
            window=wi.widget() if wi is not None else None))

    @DispatchHandle.register_listener("on_added_linelist")
    def add_linelist(self, linelist):

        # This is setting all markers at a fixed heigth in the
        # initial (before any zoom) data coordinates. Still TBD
        # how to do this in the generic case. Maybe derive heights
        # from curve data instead? Make the markers follow the
        # curve ups and downs?
        #
        # Ideally we would like to have the marker's X coordinate
        # pinned down to the plot surface in data value, and the Y
        # coordinate pinned down in screen value. This would make
        # the markers to stay at the same height in the window even
        # when the plot is zoomed. This kind of functionality doesn't
        # seem to be possible under pyqtgraph though. This requires
        # more investigation.

        plot_item = self.current_sub_window._plot_item

        # curve = plot_item.curves[0]

        data_range = plot_item.vb.viewRange()
        ymin = data_range[1][0]
        ymax = data_range[1][1]
        height = (ymax - ymin) * 0.75 + ymin

        # column names are defined in the YAML files.
        wave_column = linelist.columns['wavelength']
        id_column = linelist.columns['id']

        for i in range(len(wave_column)):
            marker = LineIDMarker(id_column[i], plot_item, orientation='vertical')

            marker.setPos(wave_column[i], height)

            plot_item.addItem(marker)
            # plot_item.addItem(marker.arrow)

            plot_item.update()

        plot_item.update()



