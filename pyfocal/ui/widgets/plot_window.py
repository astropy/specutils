from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import logging

from ..widgets.plots.plot import Plot
from .axes import DynamicAxisItem
from ...third_party.qtpy.QtWidgets import *
from ...third_party.qtpy.QtGui import *
from ..widgets.dialogs import TopAxisDialog, UnitChangeDialog

from astropy.units import Unit


class PlotWindow(QMainWindow):
    def __init__(self, **kwargs):
        super(PlotWindow, self).__init__(**kwargs)
        self._sub_window = None

        self._containers = []
        self._tool_bar = None
        self._top_axis_dialog = TopAxisDialog()
        self._unit_change_dialog = UnitChangeDialog()
        self._dynamic_axis = None
        self._plot_widget = None
        self._plot_item = None
        self._plots_units = None

    def initialize(self):
        self._dynamic_axis = DynamicAxisItem(orientation='top')
        self._plot_widget = Plot(parent=self, axisItems={'top': self._dynamic_axis})
        self.setCentralWidget(self._plot_widget)

        self._plot_item = self._plot_widget._plot_item
        self._plot_item.showAxis('top', True)
        # Add grids to the plot
        self._plot_item.showGrid(True, True)

        self._setup_connections()

    def _setup_connections(self):
        # Setup ROI connection
        act_insert_roi = self.action("actionInsert_ROI")
        act_insert_roi.triggered.connect(self._plot_widget.add_roi)

        # On accept, change the displayed axis
        self._top_axis_dialog.accepted.connect(lambda:
            self._dynamic_axis.update_axis(
                self._containers[0].layer,
                self._top_axis_dialog.ui_axis_dialog.axisModeComboBox
                    .currentIndex(),
                redshift=self._top_axis_dialog.redshift,
                ref_wave=self._top_axis_dialog.ref_wave))

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

    def add_container(self, container):
        if len(self._containers) == 0:
            self.change_units(container.layer.units[0],
                              container.layer.units[1])
        else:
            container.change_units(*self._plot_units)

        self._containers.append(container)
        self._plot_item.addItem(container.plot)

        if container.error is not None:
            self._plot_item.addItem(container.error)

        self.set_active_plot(container.layer)

        # Make sure the dynamic axis object has access to a layer
        self._dynamic_axis._layer = self._containers[0].layer

    def remove_container(self, layer):
        for container in [x for x in self._containers]:
            if container.layer == layer:
                self._plot_item.removeItem(container.plot)

                if container.error is not None:
                    self._plot_item.removeItem(container.error)

                self._containers.remove(container)

    def get_container(self, layer):
        for container in self._containers:
            if container.layer == layer:
                return container

    def change_units(self, x=None, y=None, z=None):
        for cntr in self._containers:
            cntr.change_units(x, y, z)

        self.set_labels(x_label=x, y_label=y)
        self._plot_item.enableAutoScale()
        self._plot_units = [x, y, z]

    def set_labels(self, x_label='', y_label=''):
        self._plot_item.setLabels(
            left="Flux [{}]".format(
                y_label or str(self._containers[0].layer.units[1])),
            bottom="Wavelength [{}]".format(
                x_label or str(self._containers[0].layer.units[0])))

    def set_active_plot(self, layer):
        for container in self._containers:
            if container.layer == layer:
                container.set_visibility(True, True)
            else:
                container.set_visibility(True, False)

    def set_visibility(self, layer, show, override=False):
        for container in self._containers:
            if container.layer == layer:
                container.set_visibility(show, show, inactive=False,
                                         override=override)

    def update_axis(self, layer=None, mode=None, **kwargs):
        self._dynamic_axis.update_axis(layer, mode, **kwargs)

    def _setup_toolbar_menus(self):
        # Window menu
        window_menu = QMenu()
        window_menu.addAction("Change Top Axis")
        window_menu.addAction("Change Units")

        icon = QIcon()
        icon.addPixmap(QPixmap(":/img/Settings-50.png"), QIcon.Normal,
                       QIcon.Off)

        window_menu_btn = QToolButton(self._sub_window.toolBar)
        window_menu_btn.setIcon(icon)
        window_menu_btn.setMenu(window_menu)
        window_menu_btn.setPopupMode(QToolButton.InstantPopup)

        self._sub_window.toolBar.addWidget(window_menu_btn)

        # Layer menu
        layer_menu = QMenu()
        layer_menu.addAction("Color")

        icon = QIcon()
        icon.addPixmap(QPixmap(":/img/Settings 3-50.png"), QIcon.Normal,
                        QIcon.Off)

        layer_menu_btn = QToolButton(self._sub_window.toolBar)
        layer_menu_btn.setIcon(icon)
        layer_menu_btn.setMenu(layer_menu)
        layer_menu_btn.setPopupMode(QToolButton.InstantPopup)

        self._sub_window.toolBar.addWidget(layer_menu_btn)

        # Connections
        window_menu.actions()[0].triggered.connect(
            self._top_axis_dialog.exec_)

        window_menu.actions()[1].triggered.connect(
            self._show_unit_change_dialog)

    def _show_unit_change_dialog(self):
        if self._unit_change_dialog.exec_():
            x_text = self._unit_change_dialog.disp_unit
            y_text = self._unit_change_dialog.flux_unit

            x_unit = y_unit = None

            try:
                x_unit = Unit(x_text) if x_text else None
            except ValueError as e:
                logging.error(e)

            try:
                y_unit = Unit(y_text) if y_text else None
            except ValueError as e:
                logging.error(e)

            self.change_units(x_unit, y_unit)

            self._plot_item.update()

