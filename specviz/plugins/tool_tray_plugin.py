from ..third_party.qtpy.QtWidgets import *
from ..third_party.qtpy.QtGui import *

from ..ui.widgets.plugin import Plugin
from ..ui.widgets.utils import ICON_PATH


class ToolTrayPlugin(Plugin):
    name = "Tools"
    location = "hidden"
    priority = 0

    _all_categories = {}

    def setup_ui(self):
        # ---
        # Selections setup
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
            name="Export",
            description='Export plot',
            icon_path=os.path.join(ICON_PATH, "Export-48.png"),
            category='Options',
            enabled=False)

        # ---
        # Setup button layouts
        buttons = sorted(self._tool_buttons, key=lambda x: x['category'])
        layout_horizontal = QHBoxLayout()
        self.layout_vertical.addLayout(layout_horizontal)

        for button in buttons:
            cat = self.get_category(button['category'])

            for i in range(10):
                done = False

                for j in range(4):
                    if cat['grid'].itemAtPosition(i, j):
                        continue

                    cat['grid'].addWidget(button['widget'], i, j)

                    if j < 4:
                        cat['grid'].setColumnStretch(j+1, 1)

                    done = True
                    break

                if done:
                    break

            layout_horizontal.addLayout(cat['layout'])

        # ---
        # Setup help information section
        self.label_info = QLabel()
        self.label_info.setStyleSheet("""
        QLabel {
            color: #31708f;
            background-color: #d9edf7;
            padding: 10px;
            border: 1px solid #bce8f1;
            border-radius: 4px;
        }""")
        self.label_info.setText("Hover over an icon to learn about the tool.")
        self.label_info.setWordWrap(True)

        # self.layout_vertical.addWidget(self.label_info)
        layout_horizontal.addStretch()

    def setup_connections(self):
        pass

    def get_category(self, name):
        if name is None:
            name = "user plugins"

        if name not in self._all_categories:
            vertical_layout = QVBoxLayout()
            grid_layout = QGridLayout()
            grid_layout.setSpacing(1)

            label = QLabel(self)
            label.setText(name.lower())
            label.setStyleSheet("""
                QLabel {
                    color: #999999;
                    font-variant: small-caps;
                    font-weight: bold;
                    font-size: 0.7em;
                    border-top: 1px solid #999999;
                }""")

            vertical_layout.addWidget(label)
            vertical_layout.addLayout(grid_layout)

            category = dict(name=name, grid=grid_layout, layout=vertical_layout)

            self._all_categories[name] = category

        return self._all_categories[name]

    def set_help_info(self, text):
        self.label_info.setText(text)
