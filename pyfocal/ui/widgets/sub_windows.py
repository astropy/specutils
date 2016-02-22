from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import logging
from functools import reduce

from ..widgets.plots.widgets import PlotWidget
from .axes import DynamicAxisItem
from ...third_party.qtpy.QtWidgets import *
from ...third_party.qtpy.QtGui import *
from ..widgets.dialogs import TopAxisDialog, UnitChangeDialog
from ...core.comms import Dispatch, DispatchHandle
from .region_items import LinearRegionItem

from astropy.units import Unit
import numpy as np
import pyqtgraph as pg

pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')
pg.setConfigOptions(antialias=False)


class PlotSubWindow(QMainWindow):
    def __init__(self, **kwargs):
        super(PlotSubWindow, self).__init__(**kwargs)
        self._sub_window = None

        self._containers = []
        self._tool_bar = None
        self._top_axis_dialog = TopAxisDialog()
        self._unit_change_dialog = UnitChangeDialog()
        self._dynamic_axis = None
        self._plot_widget = None
        self._plot_item = None
        self._plots_units = None
        self._rois = []
        self._equiv_width_rois = []

        DispatchHandle.setup(self)

    def initialize(self):
        self._dynamic_axis = DynamicAxisItem(orientation='top')
        self._plot_widget = pg.PlotWidget(axisItems={'top':
                                                         self._dynamic_axis})
        self.setCentralWidget(self._plot_widget)

        self._plot_item = self._plot_widget.getPlotItem()
        self._plot_item.showAxis('top', True)
        # Add grids to the plot
        self._plot_item.showGrid(True, True)

        self._setup_connections()

    def _setup_connections(self):
        # Setup ROI connection
        act_insert_roi = self.tool_bar.actions()[0]
        act_insert_roi.triggered.connect(self.add_roi)

        # On accept, change the displayed axis
        self._top_axis_dialog.accepted.connect(lambda:
            self.update_axis(
                self._containers[0].layer,
                self._top_axis_dialog.ui_axis_dialog.axisModeComboBox
                    .currentIndex(),
                redshift=self._top_axis_dialog.redshift,
                ref_wave=self._top_axis_dialog.ref_wave))

        # Setup equivalent width toggle
        act_equiv_width_mode = self.tool_bar.actions()[2]
        act_equiv_width_mode.triggered.connect(self._toggle_equiv_width)

    def _toggle_equiv_width(self, on):
        if on:
            self.add_equiv_width_rois()

            # Disable the ability to add new ROIs
            act_insert_roi = self.tool_bar.actions()[0]
            act_insert_roi.setDisabled(True)
        else:
            self.remove_equiv_width_rois()

            # Disable the ability to add new ROIs
            act_insert_roi = self.tool_bar.actions()[0]
            act_insert_roi.setDisabled(False)

    def set_sub_window(self, sub_window):
        self._sub_window = sub_window
        self._setup_toolbar_menus()

    @property
    def tool_bar(self):
        if self._tool_bar is None:
            self._tool_bar = self.findChild(QToolBar)

        return self._tool_bar

    def get_roi_mask(self, layer=None, container=None, roi=None):
        if layer is not None:
            container = self.get_container(layer)

        if container is None:
            return

        mask_holder = []
        rois = [roi] if roi is not None else self._rois

        for roi in rois:
            # roi_shape = roi.parentBounds()
            # x1, y1, x2, y2 = roi_shape.getCoords()
            x1, x2 = roi.getRegion()

            mask_holder.append((container.dispersion.value >= x1) &
                               (container.dispersion.value <= x2))

        if len(mask_holder) == 0:
            mask_holder.append(np.ones(
                shape=container.dispersion.value.shape,
                dtype=bool))

        # mask = np.logical_not(reduce(np.logical_or, mask_holder))
        mask = reduce(np.logical_or, mask_holder)
        return mask

    def add_roi(self):
        view_range = self._plot_item.viewRange()
        x_len = (view_range[0][1] - view_range[0][0]) * 0.5
        x_pos = x_len * 0.5 + view_range[0][0]

        def remove():
            self._plot_item.removeItem(roi)
            self._rois.remove(roi)

        roi = LinearRegionItem(values=[x_pos, x_pos + x_len])
        self._rois.append(roi)
        self._plot_item.addItem(roi)

        # Connect the remove functionality
        roi.sigRemoveRequested.connect(remove)

        # Connect events
        Dispatch.on_update_roi.emit(roi=roi)
        roi.sigRemoveRequested.connect(
            lambda: Dispatch.on_update_roi.emit(roi=roi))
        roi.sigRegionChangeFinished.connect(
            lambda: Dispatch.on_update_roi.emit(roi=roi))

    def add_equiv_width_rois(self):
        # First, remove existing rois
        for roi in self._rois:
            self._plot_item.removeItem(roi)

        if len(self._equiv_width_rois) == 0:
            for i in range(3):
                view_range = self._plot_item.viewRange()
                x_len = (view_range[0][1] - view_range[0][0]) * 0.25
                x_pos = view_range[0][0] + x_len * i * 1.1

                roi = LinearRegionItem(values=[x_pos, x_pos + x_len],
                                       brush=pg.mkBrush(
                                           QColor(152, 251, 152, 50)))

                if i == 1:
                    roi.setBrush(pg.mkBrush(QColor(255, 69, 0, 50)))

                self._equiv_width_rois.append(roi)

            for roi in self._equiv_width_rois:
                roi.sigRemoveRequested.connect(
                    lambda: Dispatch.on_update_roi.emit(
                        measured_rois=self._equiv_width_rois))
                roi.sigRegionChangeFinished.connect(
                    lambda: Dispatch.on_update_roi.emit(
                        measured_rois=self._equiv_width_rois))

            # Connect events
            Dispatch.on_update_roi.emit(measured_rois=self._equiv_width_rois)

        for roi in self._equiv_width_rois:
            self._plot_item.addItem(roi)

    def remove_equiv_width_rois(self):
        for roi in self._equiv_width_rois:
            self._plot_item.removeItem(roi)

        # Replace rois we removed
        for roi in self._rois:
            self._plot_item.addItem(roi)

    @DispatchHandle.register_listener("on_add_plot")
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

    @DispatchHandle.register_listener("on_remove_plot", "on_remove_layer")
    def remove_container(self, layer):
        for container in [x for x in self._containers]:
            if container.layer == layer:
                self._plot_item.removeItem(container.plot)

                if container.error is not None:
                    self._plot_item.removeItem(container.error)

                self._containers.remove(container)

    @DispatchHandle.register_listener("on_select_plot")
    def set_active_plot(self, layer):
        for container in self._containers:
            if container.layer == layer:
                container.set_visibility(True, True, inactive=False)
            else:
                container.set_visibility(True, False, inactive=True)

    def get_container(self, layer):
        for container in self._containers:
            if container.layer == layer:
                return container

    def change_units(self, x=None, y=None, z=None):
        for cntr in self._containers:
            cntr.change_units(x, y, z)

        self.set_labels(x_label=x, y_label=y)
        self._plot_item.enableAutoRange()
        self._plot_units = [x, y, z]

    def set_labels(self, x_label='', y_label=''):
        self._plot_item.setLabels(
            left="Flux [{}]".format(
                y_label or str(self._containers[0].layer.units[1])),
            bottom="Wavelength [{}]".format(
                x_label or str(self._containers[0].layer.units[0])))

    def set_visibility(self, layer, show, override=False):
        for container in self._containers:
            if container.layer == layer:
                container.set_visibility(show, show, inactive=False,
                                         override=override)

    def update_axis(self, layer=None, mode=None, **kwargs):
        self._dynamic_axis.update_axis(layer, mode, **kwargs)
        self._plot_widget.update()

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
        # layer_menu.addAction("Color")

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

    def closeEvent(self, event):
        DispatchHandle.tear_down(self)
        super(PlotSubWindow, self).closeEvent(event)

