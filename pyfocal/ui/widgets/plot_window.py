from ..widgets.plots.plot import Plot
from .axes import DynamicAxisItem

from qtpy.QtWidgets import *
from qtpy.QtGui import *

import pyqtgraph as pg


class PlotWindow(QMainWindow):
    inactive_color = pg.mkPen(color=(0, 0, 0, 75))
    active_color = pg.mkPen(color=(0, 0, 0, 255))

    def __init__(self, **kwargs):
        super(PlotWindow, self).__init__(**kwargs)
        self._sub_window = None
        self._plot_widget = None
        self._plot_item = None
        self._dynamic_axis = None

        self._containers = []
        self._tool_bar = None

    def set_sub_window(self, sub_window):
        self._sub_window = sub_window
        self._setup_toolbar_menus()

    @property
    def tool_bar(self):
        if self._tool_bar is None:
            self._tool_bar = self.findChild(QToolBar)

        return self._tool_bar

    def get_roi_mask(self, layer=None, container=None):
        if layer is not None or container is not None:
            return self._plot_widget.get_roi_mask(
                container or self.get_container(layer))

    def get_roi_data(self, layer=None, container=None):
        mask = self.get_roi_mask(layer, container)
        raise NotImplemented()

    def action(self, name):
        # TODO: Revisit this sometime in the future.
        for act in self.findChildren(QAction):
            if act.objectName() == name:
                return act

    def initialize(self):
        self._dynamic_axis = DynamicAxisItem(orientation='top')
        self._plot_widget = Plot(parent=self, axisItems={'top': self._dynamic_axis})
        self._plot_item = self._plot_widget._plot_item

        self._plot_item.showAxis('top', True)

        # Add grids to the plot
        self._plot_item.showGrid(True, True)

        self.setCentralWidget(self._plot_widget)

    def add_container(self, container):
        self._containers.append(container)

        self._plot_item.addItem(container.plot)

        if container.error is not None:
            self._plot_item.addItem(container.error)

        self.set_labels()
        self.set_active_plot(container.layer)

    def get_container(self, layer):
        for container in self._containers:
            if container.layer == layer:
                return container

    def change_unit(self, new_unit):
        for plot_container in self._containers:
            plot_container.change_unit(new_unit)

    def set_labels(self):
        self._plot_item.setLabels(
            left="Flux [{}]".format(str(self._containers[0].layer.units[1])),
            bottom="Wavelength [{}]".format(str(self._containers[
                                                    0].layer.units[0])),
        )

    def set_active_plot(self, layer):
        for container in self._containers:
            if container.layer == layer:
                container.pen = self.active_color
                container.error_pen = pg.mkPen(color=(0, 0, 0, 50))
                self.update_axis(layer=layer)
            else:
                container.pen = self.inactive_color
                container.error_pen = None

    def update_axis(self, layer=None, mode='', **kwargs):
        self._dynamic_axis.update_axis(layer, mode, **kwargs)

    def _setup_toolbar_menus(self):
        # Window menu
        window_menu = QMenu()
        window_menu.addAction("Change Top Axis")
        window_menu.addAction("Change Units")

        icon = QIcon()
        icon.addPixmap(QPixmap(":/Settings-50.png"), QIcon.Normal, QIcon.Off)

        window_menu_btn = QToolButton(self._sub_window.toolBar)
        window_menu_btn.setIcon(icon)
        window_menu_btn.setMenu(window_menu)
        window_menu_btn.setPopupMode(QToolButton.InstantPopup)

        self._sub_window.toolBar.addWidget(window_menu_btn)

        # Layer menu
        layer_menu = QMenu()
        layer_menu.addAction("Color")

        icon = QIcon()
        icon.addPixmap(QPixmap(":/Settings 3-50.png"), QIcon.Normal,
                        QIcon.Off)

        layer_menu_btn = QToolButton(self._sub_window.toolBar)
        layer_menu_btn.setIcon(icon)
        layer_menu_btn.setMenu(window_menu)
        layer_menu_btn.setPopupMode(QToolButton.InstantPopup)

        self._sub_window.toolBar.addWidget(layer_menu_btn)

