from ...third_party.qtpy.QtWidgets import *
from ...third_party.qtpy.QtGui import *


class PlotToolBar(QToolBar):
    """
    This is a separate tool bar due to the fact that you cannot add custom
    widget objects (e.g. menus) to the native tool bar within Qt Creator.
    """
    def __init__(self, *args):
        super(PlotToolBar, self).__init__(*args)
        # Disable moving
        self.setMovable(False)
        
        # Window menu
        self.window_menu = QMenu()
        self.atn_change_top_axis = self.window_menu.addAction(
            "Change Top Axis")
        self.atn_change_units = self.window_menu.addAction("Change Units")
        self.atn_line_ids = self.window_menu.addAction("Plot Line IDs")

        icon = QIcon(QPixmap(":/img/Settings-50.png"))
        self.window_menu_btn = QToolButton(self)
        self.window_menu_btn.setIcon(icon)
        self.window_menu_btn.setMenu(self.window_menu)
        self.window_menu_btn.setPopupMode(QToolButton.InstantPopup)

        self.addWidget(self.window_menu_btn)

        # Layer menu
        self.layer_menu = QMenu()
        # self.layer_menu.addAction("Color")

        icon = QIcon(QPixmap(":/img/Settings 3-50.png"))
        self.layer_menu_btn = QToolButton(self)
        self.layer_menu_btn.setIcon(icon)
        self.layer_menu_btn.setMenu(self.layer_menu)
        self.layer_menu_btn.setPopupMode(QToolButton.InstantPopup)

        self.addWidget(self.layer_menu_btn)