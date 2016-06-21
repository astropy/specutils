from ..ui.widgets.plugin import Plugin
from ..third_party.qtpy.QtWidgets import *
from ..third_party.qtpy.QtCore import *
from ..third_party.qtpy.QtGui import *

from ..ui.widgets.utils import ICON_PATH

import logging


class ToolTrayPlugin(Plugin):
    name = "Tools"
    location = "hidden"
    priority = 0

    _all_categories = {}

    def setup_ui(self):
        # ---
        # Selections setup
        self.add_tool_bar_actions(
            name="Line Labels",
            description='Add line labels',
            icon_path=os.path.join(ICON_PATH, "Label-48.png"),
            category='Selections',
            enabled=False)

        self.add_tool_bar_actions(
            name="Box ROI",
            description='Add box ROI',
            icon_path=os.path.join(ICON_PATH, "Rectangle Stroked-50.png"),
            category='Selections',
            enabled=False)

        # ---
        # Setup interactions buttons
        self.add_tool_bar_actions(
            name="Measure",
            description='Measure tool',
            icon_path=os.path.join(ICON_PATH, "Ruler-48.png"),
            category='Interactions',
            enabled=False)

        self.add_tool_bar_actions(
            name="Average",
            description='Average tool',
            icon_path=os.path.join(ICON_PATH, "Average Value-48.png"),
            category='Interactions',
            enabled=False)

        self.add_tool_bar_actions(
            name="Slice",
            description='Slice tool',
            icon_path=os.path.join(ICON_PATH, "Split Horizontal-48.png"),
            category='Interactions',
            enabled=False)

        self.add_tool_bar_actions(
            name="Detrend",
            description='Detrend tool',
            icon_path=os.path.join(ICON_PATH, "Line Chart-48.png"),
            category='Interactions',
            enabled=False)

        # ---
        # Setup transformations buttons
        self.add_tool_bar_actions(
            name="Log Scale",
            description='Log scale plot',
            icon_path=os.path.join(ICON_PATH, "Combo Chart-48.png"),
            category='Transformations',
            enabled=False)

        # ---
        # Setup plot options
        self.add_tool_bar_actions(
            name="Eport",
            description='Export plot',
            icon_path=os.path.join(ICON_PATH, "Export-48.png"),
            category='Options',
            enabled=False)

    def setup_connections(self):
        pass
